"""Do the members of one interface cluster have a similar interface? (sbatch only)

An interface cluster is a pair of chain clusters (cluster1, cluster2) of seq_cluster40_PDB, the unit MiniWorld samples PDB
interfaces by. For three kinds of interface -- protein-protein, protein-nucleic acid, antibody-antigen -- this samples interface
clusters, draws up to ``--max-members`` entries per cluster (one interface per entry, 5 A heavy-atom contact), and compares every
pair of members the way CAPRI compares a model with a reference:

* the chains of the two members are aligned by sequence (global, BLOSUM62 for protein); chain a of one member is matched to chain a
  of the other (both orders for a homo pair, the better one counts);
* each residue is one point, its pseudo-beta atom (CB, CA for Gly; C4 for A / G, C2 for C / T / U);
* reference contacts: residue pairs of the two chains with their points closer than 8 A;
* fnat      share of the reference contacts present in the other member (an unaligned residue loses its contacts);
* iRMSD     RMSD of the reference interface residues after superposing them on their counterparts;
* LRMSD     superpose the larger chain, RMSD of the smaller one;
* class     CAPRI: high (fnat >= .5 and LRMSD <= 1 or iRMSD <= 1), medium (fnat >= .3 and LRMSD <= 5 or iRMSD <= 2),
            acceptable (fnat >= .1 and LRMSD <= 10 or iRMSD <= 4), else incorrect; "unscored" when fewer than half of the
            reference interface residues are aligned (the members share a fold only in name).
"Similar" in the summary is acceptable or better. Identity is the aligned-column identity of the two chains, averaged.

usage: interface_cluster_similarity.py --interfaces DIR --db CIF_ATTACHED --out DIR [--n-clusters 3000] [--max-members 8] [--workers 56]
  DIR holds interfaces.*.tsv (scripts/sampling/extract_items.py: entry asm chain1 chain2 cluster1 cluster2 in5 seq1 seq2).
"""

import argparse
import glob
import json
import multiprocessing as mp
import os
import random
import sys
from collections import defaultdict
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"

import lmdb  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402

sys.path.insert(0, str(Path(__file__).parents[1] / "design"))
from build_stats import shard_paths  # noqa: E402

from structcooker.instructions.mols import CIFMolAttached  # noqa: E402
from structcooker.instructions.potential.pairs import (  # noqa: E402
    AA,
    TYPES,
    representative_atoms,
    residue_types,
)
from structcooker.instructions.transforms.codecs import from_bytes  # noqa: E402

KINDS = ("protein_protein", "protein_na", "antibody_antigen")
CONTACT = 8.0
_ST: dict = {}


def kind_of(k1, k2):
    t = {k1[5], k2[5]}
    if t == {"P"}:
        return "protein_protein"
    if len(t) == 2 and "P" in t and t <= {"P", "R", "D"}:
        return "protein_na"
    if t == {"A", "P"}:
        return "antibody_antigen"
    return None


def sample(args):
    groups = defaultdict(dict)  # (k1, k2) -> {entry: (asm, chain of k1, chain of k2)}
    for p in sorted(glob.glob(str(Path(args.interfaces) / "interfaces.*.tsv"))):
        with open(p) as f:
            for ln in f:
                e, asm, c1, c2, k1, k2, in5 = ln.split("\t", 7)[:7]
                if in5 != "1":
                    continue
                if k1 > k2:
                    k1, k2, c1, c2 = k2, k1, c2, c1
                if kind_of(k1, k2) is None:
                    continue
                groups[(k1, k2)].setdefault(e, (asm, c1, c2))
    rng = random.Random(0)
    by_kind = defaultdict(list)
    for key, ents in groups.items():
        if len(ents) >= 2:
            by_kind[kind_of(*key)].append(key)
    tasks, stats = [], {}
    for kind in KINDS:
        keys = sorted(by_kind[kind])
        stats[kind] = {
            "interface_clusters": sum(
                1 for k, e in groups.items() if kind_of(*k) == kind
            ),
            "with_2plus_entries": len(keys),
            "clusters_above_max_members": sum(
                len(groups[k]) > args.max_members for k in keys
            ),
        }
        rng.shuffle(keys)
        for key in keys[: args.n_clusters]:
            ents = sorted(groups[key].items())
            rng.shuffle(ents)
            members = [(e, *v) for e, v in ents[: args.max_members]]
            tasks.append((kind, key, len(groups[key]), members))
    print(json.dumps(stats), flush=True)
    return tasks, stats


def _init(db):
    _ST["envs"] = [
        lmdb.open(
            p,
            readonly=True,
            lock=False,
            readahead=False,
            subdir=Path(p).is_dir(),
            max_readers=512,
        )
        for p in shard_paths(db)
    ]
    _ST["aln"] = {}
    for name, mat in (("prot", substitution_matrices.load("BLOSUM62")), ("na", None)):
        a = Align.PairwiseAligner()
        a.mode = "global"
        if mat is not None:
            a.substitution_matrix, a.open_gap_score, a.extend_gap_score = (
                mat,
                -10.0,
                -0.5,
            )
        else:
            a.match_score, a.mismatch_score, a.open_gap_score, a.extend_gap_score = (
                2.0,
                -1.0,
                -5.0,
                -1.0,
            )
        _ST["aln"][name] = a


def load_chain(cm, rtype, rep, chain_of_res, name):
    ids = cm.chains.chain_id.value
    ci = np.flatnonzero(ids == name)
    if len(ci) == 0:
        return None
    res = np.flatnonzero((chain_of_res == ci[0]) & (rtype >= 0))
    if len(res) == 0:
        return None
    xyz = np.full((len(res), 3), np.nan)
    ok = rep[res] >= 0
    xyz[ok] = cm.atoms.xyz.value[rep[res][ok]]
    na = bool((rtype[res] > len(AA)).any())
    seq = "".join(TYPES[t][1] if t > len(AA) else TYPES[t] for t in rtype[res].tolist())
    return {"xyz": xyz, "seq": seq, "na": na}


def load_member(entry, asm, ca, cb):
    raw = None
    for env in _ST["envs"]:
        with env.begin() as t:
            raw = t.get(entry.encode())
        if raw is not None:
            break
    cm = CIFMolAttached.from_dict(from_bytes(bytes(raw))[asm]["cifmol_attached_dict"])
    rtype, _ = residue_types(cm)
    rep = representative_atoms(cm, rtype)
    chain_of_res = cm.index_table.residues_to_chains(np.arange(len(rtype)))
    a, b = (
        load_chain(cm, rtype, rep, chain_of_res, ca),
        load_chain(cm, rtype, rep, chain_of_res, cb),
    )
    return None if a is None or b is None else (a, b)


def align(x, y):
    """(index in x, index in y) of the aligned columns, and their identity."""
    aln = _ST["aln"]["na" if x["na"] and y["na"] else "prot"]
    best = aln.align(x["seq"], y["seq"])[0]
    ix, iy = [], []
    for (a0, a1), (b0, b1) in zip(*best.aligned, strict=True):
        ix.extend(range(a0, a1))
        iy.extend(range(b0, b1))
    ix, iy = np.array(ix, dtype=int), np.array(iy, dtype=int)
    ident = (
        float(
            np.mean([x["seq"][i] == y["seq"][j] for i, j in zip(ix, iy, strict=True)])
        )
        if len(ix)
        else 0.0
    )
    return ix, iy, ident


def cdist(a, b):
    return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))


def kabsch(p, q):
    """Rotation, translation fitting p onto q (rows are points); returns the fitted p."""
    pc, qc = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    return (p - pc) @ rot + qc, (pc, qc, rot)


def rmsd(p, q):
    return float(np.sqrt(((p - q) ** 2).sum(1).mean()))


def compare(ref, mod):
    """CAPRI measures of ``mod`` against ``ref`` (each ``(chain a, chain b)``), chains already matched a-a, b-b."""
    ra, rb = ref
    ma, mb = mod
    ixa, iya, ida = align(ra, ma)
    ixb, iyb, idb = align(rb, mb)
    mapa = np.full(len(ra["seq"]), -1)
    mapb = np.full(len(rb["seq"]), -1)
    mapa[ixa], mapb[ixb] = iya, iyb
    d = cdist(ra["xyz"], rb["xyz"])
    contact = d < CONTACT
    n_contact = int(contact.sum())
    if n_contact == 0:
        return None
    ia, ib = np.flatnonzero(contact.any(1)), np.flatnonzero(contact.any(0))
    ok_a = ia[(mapa[ia] >= 0)]
    ok_b = ib[(mapb[ib] >= 0)]
    ok_a = ok_a[np.isfinite(ma["xyz"][mapa[ok_a]]).all(1)]
    ok_b = ok_b[np.isfinite(mb["xyz"][mapb[ok_b]]).all(1)]
    cov = (len(ok_a) + len(ok_b)) / (len(ia) + len(ib))
    out = {
        "cov": cov,
        "identity": (ida + idb) / 2,
        "fnat": np.nan,
        "irmsd": np.nan,
        "lrmsd": np.nan,
        "cls": "unscored",
    }
    if cov < 0.5 or len(ok_a) + len(ok_b) < 3:
        return out
    pi = np.vstack([ra["xyz"][ok_a], rb["xyz"][ok_b]])
    qi = np.vstack([ma["xyz"][mapa[ok_a]], mb["xyz"][mapb[ok_b]]])
    out["irmsd"] = rmsd(kabsch(pi, qi)[0], qi)
    pa, pb = np.flatnonzero(mapa >= 0), np.flatnonzero(mapb >= 0)
    pa = pa[np.isfinite(ra["xyz"][pa]).all(1) & np.isfinite(ma["xyz"][mapa[pa]]).all(1)]
    pb = pb[np.isfinite(rb["xyz"][pb]).all(1) & np.isfinite(mb["xyz"][mapb[pb]]).all(1)]
    dm = cdist(ma["xyz"], mb["xyz"])
    ca = contact[np.ix_(pa, pb)] if len(pa) and len(pb) else np.zeros((0, 0), bool)
    kept = 0
    if ca.size:
        sub = dm[np.ix_(mapa[pa], mapb[pb])] < CONTACT
        kept = int((ca & sub).sum())
    out["fnat"] = kept / n_contact
    big, small = (0, 1) if len(ra["seq"]) >= len(rb["seq"]) else (1, 0)
    sel = {0: (pa, mapa, ra, ma), 1: (pb, mapb, rb, mb)}
    ib_, mp_, rr, mm = sel[big]
    isl, ms_, rs_, mms = sel[small]
    if len(ib_) >= 3 and len(isl) >= 1:
        _, (pc, qc, rot) = kabsch(rr["xyz"][ib_], mm["xyz"][mp_[ib_]])
        moved = (rs_["xyz"][isl] - pc) @ rot + qc
        out["lrmsd"] = rmsd(moved, mms["xyz"][ms_[isl]])
    f, i, l = out["fnat"], out["irmsd"], out["lrmsd"]
    lo = lambda x, t: bool(np.isfinite(x) and x <= t)
    if f >= 0.5 and (lo(l, 1) or lo(i, 1)):
        out["cls"] = "high"
    elif f >= 0.3 and (lo(l, 5) or lo(i, 2)):
        out["cls"] = "medium"
    elif f >= 0.1 and (lo(l, 10) or lo(i, 4)):
        out["cls"] = "acceptable"
    else:
        out["cls"] = "incorrect"
    return out


RANK = {"high": 3, "medium": 2, "acceptable": 1, "incorrect": 0, "unscored": -1}


def work(task):
    kind, (k1, k2), n_entries, members = task
    loaded = []
    for e, asm, ca, cb in members:
        try:
            m = load_member(e, asm, ca, cb)
        except Exception:  # noqa: BLE001  (a unit that cannot be loaded is skipped)
            m = None
        loaded.append((e, m))
    rows = []
    for i in range(len(loaded)):
        for j in range(i + 1, len(loaded)):
            (ei, mi), (ej, mj) = loaded[i], loaded[j]
            if mi is None or mj is None:
                continue
            options = [compare(mi, mj)]
            if k1 == k2:
                options.append(compare(mi, (mj[1], mj[0])))
            options = [o for o in options if o is not None]
            if not options:
                continue
            best = max(
                options,
                key=lambda o: (
                    RANK[o["cls"]],
                    -(o["irmsd"] if np.isfinite(o["irmsd"]) else 99),
                ),
            )
            rows.append(
                (
                    kind,
                    k1,
                    k2,
                    n_entries,
                    ei,
                    ej,
                    best["identity"],
                    best["cov"],
                    best["fnat"],
                    best["irmsd"],
                    best["lrmsd"],
                    best["cls"],
                )
            )
    return rows


COLS = [
    "kind",
    "cluster1",
    "cluster2",
    "n_entries",
    "entry_i",
    "entry_j",
    "identity",
    "coverage",
    "fnat",
    "irmsd",
    "lrmsd",
    "class",
]


def summarize(path, out):
    import csv

    rows = list(csv.DictReader(open(path), delimiter="\t"))
    for r in rows:
        r["identity"], r["n_entries"] = float(r["identity"]), int(r["n_entries"])
        r["id_bin"] = (
            "<60%"
            if r["identity"] < 0.6
            else "60-90%"
            if r["identity"] < 0.9
            else ">=90%"
        )
        n = r["n_entries"]
        r["size_bin"] = (
            "2" if n == 2 else "3-5" if n <= 5 else "6-20" if n <= 20 else ">20"
        )
    lines = [
        "# Are members of one interface cluster similar?",
        "",
        "CAPRI-style comparison of two members of the same interface cluster (pseudo-beta points, 8 A contacts; see the script header).",
        "",
    ]

    def table(title, by, order):
        lines.extend(
            [
                f"## {title}",
                "",
                f"| {by or 'kind'} | pairs | clusters | unscored | incorrect | acceptable | medium | high | similar (acceptable or better) | similar, of scored pairs | medium or better |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for kind in KINDS:
            for lab in order:
                g = [
                    r
                    for r in rows
                    if r["kind"] == kind and (by is None or r[by] == lab)
                ]
                if not g:
                    continue
                n = len(g)
                c = {k: sum(r["class"] == k for r in g) for k in RANK}
                sc = n - c["unscored"]
                sim = c["acceptable"] + c["medium"] + c["high"]
                name = kind if by is None else f"{kind} / {lab}"
                lines.append(
                    f"| {name} | {n:,} | {len({(r['cluster1'], r['cluster2']) for r in g}):,} | {c['unscored'] / n:.1%} | {c['incorrect'] / n:.1%} | {c['acceptable'] / n:.1%} | "
                    f"{c['medium'] / n:.1%} | {c['high'] / n:.1%} | {sim / n:.1%} | {(sim / sc if sc else float('nan')):.1%} | {(c['medium'] + c['high']) / n:.1%} |"
                )
        lines.append("")

    table(
        "All pairs (each interface cluster contributes up to C(members, 2) pairs)",
        None,
        [None],
    )
    table(
        "By sequence identity of the aligned chains",
        "id_bin",
        ["<60%", "60-90%", ">=90%"],
    )
    table(
        "By number of entries in the interface cluster",
        "size_bin",
        ["2", "3-5", "6-20", ">20"],
    )
    lines.extend(
        [
            "## Per cluster (each cluster counts once)",
            "",
            "| kind | clusters | mean share of similar pairs | clusters with every pair similar | clusters with no pair similar |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for kind in KINDS:
        per = {}
        for r in rows:
            if r["kind"] == kind:
                per.setdefault((r["cluster1"], r["cluster2"]), []).append(
                    r["class"] in ("acceptable", "medium", "high")
                )
        if per:
            rate = [sum(v) / len(v) for v in per.values()]
            lines.append(
                f"| {kind} | {len(per):,} | {np.mean(rate):.1%} | {np.mean([x == 1 for x in rate]):.1%} | {np.mean([x == 0 for x in rate]):.1%} |"
            )
    Path(out).write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interfaces", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-clusters", type=int, default=3000)
    ap.add_argument("--max-members", type=int, default=8)
    ap.add_argument("--workers", type=int, default=56)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tasks, stats = sample(a)
    (out / "sampling.json").write_text(json.dumps(stats, indent=1))
    n = 0
    with (
        open(out / "pairs.tsv", "w") as f,
        mp.get_context("fork").Pool(
            a.workers, initializer=_init, initargs=(a.db,)
        ) as pool,
    ):
        f.write("\t".join(COLS) + "\n")
        for k, rows in enumerate(pool.imap_unordered(work, tasks, chunksize=4), 1):
            for r in rows:
                f.write(
                    "\t".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in r)
                    + "\n"
                )
                n += 1
            if k % 500 == 0:
                print(f"{k}/{len(tasks)} clusters, {n} pairs", flush=True)
    print(f"done: {len(tasks)} clusters, {n} pairs", flush=True)
    summarize(out / "pairs.tsv", out / "summary.md")


if __name__ == "__main__":
    main()
