"""Find the PDB entries and chains that carry a cyclic peptide, a PTM residue or a glycan (sbatch for a full scan).

usage: pdb_special_chains.py --out DIR [--ids 1cwa,1nb1,...] [--db cif_pdb.lmdb] [--workers N] [--every K]
                             [--max-pep-len 50] [--ccd-only]

Reads `lmdb/pdb/cif/cif_pdb.lmdb` (schema A, one record per PDB id; the first assembly, first model of each) and the raw CCD
(`materials/raw/ccd/components_*.cif.gz`, whose `_chem_comp.type` / `mon_nstd_parent_comp_id` / `pdbx_pcm` the CCD LMDB does not
keep). Writes into DIR:
  ccd_classes.tsv      comp_id, type, parent, pcm, class  (class: ptm_L | ncaa_D | sugar | other)
  special_chains.tsv   pdb_id, chain (auth_asym_id), category, chain_len, evidence
  entries.tsv          pdb_id, assembly key, categories
  summary.txt          counts per category (chains, entries)

Categories (the notions follow RFD4-Proteina's PDB components):
  ptm                  a polypeptide chain with a residue whose CCD component is peptide-linking and has a canonical amino
                       acid as `mon_nstd_parent_comp_id` (RFD4 `pdb_ptm`; MSE excluded). D-amino acids (parent L form) are
                       reported apart as ncaa_d_parent. evidence = comp:count list, `pcm` marks the CCD's own PTM flag.
  glycan_branched      an entity of type `branched` (RFD4 `pdb_glycan`, branched components); evidence = sugars, attachment
  glycan_single        a non-polymer chain whose CCD type is a saccharide (RFD4: standalone saccharides)
  glycosylated_chain   a polypeptide chain with a residue covalently bonded (struct_conn covale) to a sugar residue
  cyclic_headtail      polypeptide with an N(first)-C(last) closure: a struct_conn covale, or, failing that, N-C <= 1.8 A
  cyclic_isopeptide    polypeptide with a side-chain N-C amide bond between residues >= 2 apart
  cyclic_disulfide     polypeptide of <= --max-pep-len residues with an intra-chain disulfide (RFD4 `cpsea_pdb` closure types)
  cyclic_other         any other intra-chain covale between residues >= 2 apart (thioether, lanthionine, ...), same length cap
A chain is reported once per category (symmetry copies share their auth_asym_id).
"""

from __future__ import annotations

import argparse
import gzip
import multiprocessing as mp
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import lmdb
import numpy as np

from structcooker.instructions.transforms.codecs import from_bytes

ROOT = Path("/data/shared/cssb_data/BioMol_clean")
DEFAULT_DB = ROOT / "lmdb/pdb/cif/cif_pdb.lmdb"
CCD_GZ = ROOT / "materials/raw/ccd/components_20260803.cif.gz"
CANONICAL = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
TIERS = (
    (20_000_000, "workers"),
    (200_000_000, "mid"),
    (None, "huge"),
)  # value bytes -> pool; decoded records are ~100x larger
HEADTAIL_MAX_A = 1.8
_S: dict = {}


# ---- CCD classes -------------------------------------------------------------------------------
def build_ccd_classes(gz: Path) -> dict[str, tuple[str, str, str, str]]:
    """comp_id -> (type, parent, pcm, class); items of each data_ block precede its first loop_."""
    out: dict[str, tuple[str, str, str, str]] = {}
    cur: dict[str, str] = {}
    open_block = False

    def flush():
        cid = cur.get("id")
        if not cid:
            return
        typ = cur.get("type", "").upper()
        parent = cur.get("mon_nstd_parent_comp_id", "?")
        parents = [
            p for p in re.split(r"[\s,]+", parent.strip("\"'")) if p and p != "?"
        ]
        if "SACCHARIDE" in typ:
            cls = "sugar"
        elif (
            "PEPTIDE" in typ
            and "PEPTIDE-LIKE" not in typ
            and any(p in CANONICAL for p in parents)
        ):
            if cid in CANONICAL or cid == "MSE":
                cls = "other"
            elif typ.startswith("D-"):
                cls = "ncaa_D"
            else:
                cls = "ptm_L"
        else:
            cls = "other"
        out[cid] = (typ, ",".join(parents) or "?", cur.get("pdbx_pcm", "?"), cls)

    with gzip.open(gz, "rt") as f:
        for line in f:
            if line.startswith("data_"):
                if cur:
                    flush()
                cur, open_block = {}, True
            elif open_block:
                if line.startswith(("loop_", "#")):
                    if line.startswith("loop_"):
                        open_block = False
                    continue
                if line.startswith("_chem_comp."):
                    k, _, v = line[len("_chem_comp.") :].partition(" ")
                    if k in ("id", "type", "mon_nstd_parent_comp_id", "pdbx_pcm"):
                        cur[k] = v.strip().strip("\"'")
    if cur:
        flush()
    return out


# ---- per-record detection ------------------------------------------------------------------------
def _pick_assembly(a: dict) -> str:
    def k(s):
        return [int(x) if x.isdigit() else 0 for x in s.split("_")]

    return sorted(a, key=k)[0]


def _fmt(counter: Counter, limit=8) -> str:
    items = counter.most_common(limit)
    s = ",".join(f"{c}:{n}" for c, n in items)
    return s + (f",+{len(counter) - limit}" if len(counter) > limit else "")


def detect(
    pdb_id: str, rec: dict, ccd_cls: dict, max_pep_len: int
) -> tuple[str, list[tuple]]:
    a = rec["assembly_dict"]
    if (
        not a
    ):  # an entry whose record holds no assembly (32 of 249,676): nothing to scan
        return "", []
    akey = _pick_assembly(a)
    g = a[akey]
    A, R, C, ix = g["atoms"], g["residues"], g["chains"], g["index_table"]
    atom_name = A["nodes"]["id"]["value"].astype(str)
    element = A["nodes"]["element"]["value"].astype(str)
    xyz = A["nodes"]["xyz"]["value"]
    ccd = R["nodes"]["chem_comp_id"]["value"].astype(str)
    auth_idx = R["nodes"]["auth_idx"]["value"]
    etype = C["nodes"]["entity_type"]["value"].astype(str)
    auth = C["nodes"]["auth_asym_id"]["value"].astype(str)
    a2r = np.asarray(ix["atom_to_res"])
    r2c = np.asarray(ix["res_to_chain"])
    cptr = np.asarray(ix["chain_res_indptr"])
    cres = np.asarray(ix["chain_res_indices"])
    rptr = np.asarray(ix["res_atom_indptr"])
    ratoms = np.asarray(ix["res_atom_indices"])
    cls_of = (
        np.array([ccd_cls.get(c, ("", "", "", "other"))[3] for c in ccd])
        if len(ccd)
        else ccd
    )
    pcm_of = np.array([ccd_cls.get(c, ("", "", "?", ""))[2] == "Y" for c in ccd])
    pos = np.zeros(len(ccd), dtype=np.int64)  # residue index -> position in its chain
    for c in range(len(etype)):
        r = cres[cptr[c] : cptr[c + 1]]
        pos[r] = np.arange(len(r))

    rows: dict[
        tuple[str, str], list
    ] = {}  # (chain, category) -> [chain_len, evidence parts]

    def add(chain_i: int, cat: str, ev: str):
        key = (auth[chain_i], cat)
        n = int(cptr[chain_i + 1] - cptr[chain_i])
        if key in rows:
            if ev and ev not in rows[key][1]:
                rows[key][1].append(ev)
        else:
            rows[key] = [n, [ev] if ev else []]

    def label(ri: int) -> str:
        return f"{auth[r2c[ri]]}:{ccd[ri]}{auth_idx[ri]}"

    # struct_conn covale / disulf pairs (atom indices -> residues)
    sc = A["edges"].get("struct_conn")
    covale, disulf = [], []
    if sc is not None and len(sc["src_indices"]):
        kinds = np.asarray(sc["value"]).astype(str).reshape(-1, 2)[:, 0]
        for i, j, k in zip(sc["src_indices"], sc["dst_indices"], kinds, strict=False):
            if k == "covale":
                covale.append((int(i), int(j)))
            elif k == "disulf":
                disulf.append((int(i), int(j)))

    # sugar attachments (covale between a sugar residue and anything else)
    attach: dict[int, list[str]] = defaultdict(
        list
    )  # sugar chain index -> "->A:ASN123"
    for i, j in covale:
        ri, rj = int(a2r[i]), int(a2r[j])
        for s, o, sa, oa in ((ri, rj, i, j), (rj, ri, j, i)):
            if cls_of[s] == "sugar" and cls_of[o] != "sugar" and r2c[s] != r2c[o]:
                attach[r2c[s]].append(f"->{label(o)}({atom_name[oa]}-{atom_name[sa]})")
                if etype[r2c[o]].startswith("polypeptide"):
                    add(
                        r2c[o],
                        "glycosylated_chain",
                        f"{ccd[o]}{auth_idx[o]}-{ccd[s]}@{auth[r2c[s]]}",
                    )

    for c in range(len(etype)):
        res = cres[cptr[c] : cptr[c + 1]]
        n = len(res)
        if etype[c] == "branched":
            sug = Counter(ccd[res].tolist())
            att = attach.get(c, [])
            add(
                c,
                "glycan_branched",
                f"{_fmt(sug)};{'|'.join(dict.fromkeys(att)) or 'unattached'}",
            )
        elif etype[c] == "non-polymer":
            if n and cls_of[res[0]] == "sugar":
                att = attach.get(c, [])
                add(
                    c,
                    "glycan_single",
                    f"{ccd[res[0]]};{'|'.join(dict.fromkeys(att)) or 'unattached'}",
                )
        elif etype[c].startswith("polypeptide"):
            for cat, name in (("ptm_L", "ptm"), ("ncaa_D", "ncaa_d_parent")):
                m = cls_of[res] == cat
                if m.any():
                    cnt = Counter(ccd[res][m].tolist())
                    npcm = int(pcm_of[res][m].sum())
                    add(c, name, f"{_fmt(cnt)};n={int(m.sum())};pcm={npcm}")

    # cyclic closures, per polypeptide chain
    pep = [c for c in range(len(etype)) if etype[c].startswith("polypeptide")]
    pep_set = set(pep)
    closure: dict[int, set] = defaultdict(set)
    for i, j in covale:
        ri, rj = int(a2r[i]), int(a2r[j])
        c = int(r2c[ri])
        if c != int(r2c[rj]) or c not in pep_set:
            continue
        L = int(cptr[c + 1] - cptr[c])
        pi, pj = int(pos[ri]), int(pos[rj])
        ni, nj = atom_name[i], atom_name[j]
        if {ni, nj} == {"N", "C"}:
            n_pos, c_pos = (pi, pj) if ni == "N" else (pj, pi)
            if n_pos == 0 and c_pos == L - 1 and L > 1:
                add(
                    c,
                    "cyclic_headtail",
                    f"N{auth_idx[ri if ni == 'N' else rj]}-C{auth_idx[rj if ni == 'N' else ri]};bond",
                )
                closure[c].add("ht")
                continue
        if abs(pi - pj) < 2:
            continue
        pair = f"{label(ri)}.{ni}-{label(rj)}.{nj}"
        if {element[i], element[j]} == {"N", "C"} and not (
            ni in ("N", "C") and nj in ("N", "C")
        ):
            add(c, "cyclic_isopeptide", pair)
        elif max_pep_len >= L:
            add(c, "cyclic_other", pair)
    for i, j in disulf:
        ri, rj = int(a2r[i]), int(a2r[j])
        c = int(r2c[ri])
        if (
            c == int(r2c[rj])
            and c in pep_set
            and int(cptr[c + 1] - cptr[c]) <= max_pep_len
            and abs(pos[ri] - pos[rj]) >= 2
        ):
            add(
                c,
                "cyclic_disulfide",
                f"{ccd[ri]}{auth_idx[ri]}-{ccd[rj]}{auth_idx[rj]}",
            )
    for c in pep:  # geometric head-to-tail when the entry has no closure bond
        L = int(cptr[c + 1] - cptr[c])
        if L < 3 or "ht" in closure[c]:
            continue
        r0, r1 = int(cres[cptr[c]]), int(cres[cptr[c + 1] - 1])
        a0 = [x for x in ratoms[rptr[r0] : rptr[r0 + 1]] if atom_name[x] == "N"]
        a1 = [x for x in ratoms[rptr[r1] : rptr[r1 + 1]] if atom_name[x] == "C"]
        if a0 and a1:
            d = float(np.linalg.norm(xyz[a0[0]] - xyz[a1[0]]))
            if np.isfinite(d) and d <= HEADTAIL_MAX_A:
                add(
                    c,
                    "cyclic_headtail",
                    f"N{auth_idx[r0]}-C{auth_idx[r1]};geom={d:.2f}A",
                )

    out = [
        (pdb_id, ch, cat, v[0], ";".join(v[1])[:400]) for (ch, cat), v in rows.items()
    ]
    return akey, out


# ---- driver ------------------------------------------------------------------------------------------
def _init(db, ccd_cls, max_pep_len):
    _S["env"] = lmdb.open(
        str(db),
        readonly=True,
        lock=False,
        readahead=False,
        max_readers=512,
        subdir=Path(db).is_dir(),
    )
    _S["ccd"], _S["max"] = ccd_cls, max_pep_len


def _work(key):
    with _S["env"].begin() as txn:
        raw = txn.get(key.encode())
    if raw is None:
        return key, None, [], "missing"
    try:
        akey, rows = detect(key, from_bytes(raw), _S["ccd"], _S["max"])
        return key, akey, rows, ""  # noqa: TRY300
    except Exception as e:  # noqa: BLE001  keep the scan going; the error is recorded per entry
        import traceback

        return (
            key,
            None,
            [],
            f"{type(e).__name__}: {e} @{traceback.extract_tb(e.__traceback__)[-1].lineno}"[
                :200
            ],
        )


def read_index(db: Path):
    idx = Path(str(db) + ".index.tsv")
    keys = []
    for ln in idx.read_text().splitlines():
        if ln.startswith(("#", "key\t")):
            continue
        k, size = ln.split("\t")[:2]
        keys.append((k, int(size)))
    return keys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--ccd", default=str(CCD_GZ))
    ap.add_argument(
        "--ids", help="comma-separated PDB ids (test); default: every key of the DB"
    )
    ap.add_argument("--workers", type=int, default=96)
    ap.add_argument("--mid-workers", type=int, default=8)
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--max-pep-len", type=int, default=50)
    ap.add_argument("--ccd-only", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tsv = out / "ccd_classes.tsv"
    if tsv.exists():
        ccd_cls = {}
        for ln in tsv.read_text().splitlines()[1:]:
            c, t, p, pcm, cls = ln.split("\t")
            ccd_cls[c] = (t, p, pcm, cls)
    else:
        ccd_cls = build_ccd_classes(Path(args.ccd))
        with open(tsv, "w") as f:
            f.write("comp_id\ttype\tparent\tpcm\tclass\n")
            for c, (t, p, pcm, cls) in sorted(ccd_cls.items()):
                f.write(f"{c}\t{t}\t{p}\t{pcm}\t{cls}\n")
    print("CCD classes:", dict(Counter(v[3] for v in ccd_cls.values())), flush=True)
    if args.ccd_only:
        return

    db = Path(args.db)
    if args.ids:
        keys = [(k, 0) for k in args.ids.split(",")]
    else:
        keys = read_index(db)[:: args.every]
    small = [k for k, s in keys if s <= TIERS[0][0]]
    mid = [k for k, s in keys if TIERS[0][0] < s <= TIERS[1][0]]
    big = [k for k, s in keys if s > TIERS[1][0]]
    print(f"{len(small)} small + {len(mid)} mid + {len(big)} huge records", flush=True)

    results = []
    initargs = (db, ccd_cls, args.max_pep_len)
    for ks, w, chunk in (
        (small, args.workers, 64),
        (mid, args.mid_workers, 1),
        (big, 1, 1),
    ):
        if not ks:
            continue
        with mp.Pool(
            w, _init, initargs, maxtasksperchild=None if chunk > 1 else 1
        ) as pool:
            for n, r in enumerate(pool.imap_unordered(_work, ks, chunksize=chunk), 1):
                results.append(r)
                if n % 20000 == 0:
                    print(f"  {n}/{len(ks)}", flush=True)

    results.sort(key=lambda r: r[0])
    cat_chains, cat_entries, errors = Counter(), defaultdict(set), []
    with (
        open(out / "special_chains.tsv", "w") as f,
        open(out / "entries.tsv", "w") as e,
    ):
        f.write("pdb_id\tchain\tcategory\tchain_len\tevidence\n")
        e.write("pdb_id\tassembly\tcategories\n")
        for key, akey, rows, err in results:
            if err:
                errors.append((key, err))
                continue
            cats = sorted({r[2] for r in rows})
            for r in sorted(rows, key=lambda r: (r[1], r[2])):
                f.write("\t".join(map(str, r)) + "\n")
                cat_chains[r[2]] += 1
                cat_entries[r[2]].add(key)
            if cats:
                e.write(f"{key}\t{akey}\t{','.join(cats)}\n")
    with open(out / "summary.txt", "w") as f:
        f.write(
            f"entries scanned: {len(results) - len(errors)} (errors: {len(errors)})\n"
        )
        f.write("category\tchains\tentries\n")
        for c in sorted(cat_chains):
            f.write(f"{c}\t{cat_chains[c]}\t{len(cat_entries[c])}\n")
        for k, m in errors[:50]:
            f.write(f"ERROR\t{k}\t{m}\n")
    print((out / "summary.txt").read_text(), flush=True)


if __name__ == "__main__":
    sys.exit(main())
