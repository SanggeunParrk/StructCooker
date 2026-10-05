"""Stage the AFDB-distributed MSAs the AFM structures need into materials/raw/afdb_multimer/msa.

The release ships one MSA per MONOMER entity (msas/251124_1.2M, 251124_8.8M, 251209_13.4M;
23,425,744 entities). A homodimer uses the MSA of its own entity id. A heterodimer uses the
MSAs of its two monomers -- the paper concatenates them without pairing -- so each chain has
to be tied to a monomer entity, which the release does not do for us:

* chain A <-> ``uniprot_ac_1``, chain B <-> ``uniprot_ac_2`` of heterodimer_metadata.csv;
* the UniProt accession -> monomer entity comes from homodimer_metadata.csv;
* 2,229 chains have no such row (their monomers are among the ~2M distributed MSAs that
  homodimer_metadata.csv never lists), so they were tied by exact sequence instead
  (``msa_found_by_sequence.tsv``, measured 2026-09-25).

Because that link is derived, every heterodimer chain's MSA is re-checked here: its query
sequence must equal the chain's own sequence in the structure. The resulting table is ours, not the release's, so it goes to
``BioMol/materials/intermediate/afdb_multimer/heterodimer/chain_msa.tsv`` (raw holds only
what was downloaded). Generated MSAs (entities the release does not ship) live in
``intermediate/afdb_multimer/msa_generated/`` as ``.a3m`` and, compressed, ``.a3m.zst``.

Only the ~1.9M MSAs actually used are staged. Directory-stored ones are HARD-LINKED (same
filesystem: no extra space, and the raw copy survives if the download mirror is removed);
tar-only ones (the 8.8M batch) are extracted member by member with GNU tar -- Python's
tarfile ends up reading every member's data on this filesystem (397 GB in 22 min).

Subcommands (run as separate batch jobs):
    plan                         build the work list and the heterodimer chain map
    stage --task T --ntasks N    link / extract one slice (array job)
    finalize                     add generated MSAs, verify, write SOURCE.tsv + README
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

_DATA = Path(os.environ.get("DATA_ROOT", "/data/shared/cssb_data"))
# The download (release msas/ batches + metadata CSVs) and a scratch dir for the plan.
DIST = Path(os.environ.get("AFM_MSA_DOWNLOAD", "/data/psk6950/external_source/AFDB/multimer/msas"))
BATCHES = ("251124_1.2M", "251124_8.8M", "251209_13.4M")
WORK = Path(os.environ.get("AFM_STAGE_WORK", "/data/psk6950/external_source/AFDB/multimer"))
RAW = _DATA / "BioMol/materials/raw/afdb_multimer"            # release MSAs, as shipped
INTER = _DATA / "BioMol/materials/intermediate/afdb_multimer"  # what we derive from them
HOMO_META = WORK / "homodimer_metadata.csv"
HET_META = Path(os.environ.get("AFM_HET_METADATA", "/data/shared/cssb_data/AFDB_heterodimer/heterodimer_metadata.csv"))
SEQ_FOUND = WORK / "msa_found_by_sequence.tsv"
PLAN = WORK / "msa_stage_plan.tsv"          # entity, kind(dir|tar), path, member, roles
CHAIN_SEQ = WORK / "msa_stage_chainseq.tsv"  # het, chain, sequence

_MSA = re.compile(r"(AF[-_]\d+)-msa_v1\.a3m\.zst$")
_AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def _norm(x: str) -> str:
    return x.replace("AF_", "AF-")


def _dest(entity: str) -> Path:
    n = entity.split("-", 1)[1]
    return RAW / "msa" / n[-3:] / f"{entity}-msa_v1.a3m.zst"


# ---------------------------------------------------------------- plan
def _inventory() -> dict[str, tuple[str, str, str]]:
    """Map entity -> (kind, path, member); directory copies win over tar copies."""
    inv: dict[str, tuple[str, str, str]] = {}
    tars: list[str] = []
    for batch in BATCHES:
        for entry in os.scandir(DIST / batch):
            if entry.is_dir():
                for f in os.scandir(entry.path):
                    m = _MSA.search(f.name)
                    if m:
                        inv.setdefault(_norm(m.group(1)), ("dir", f.path, ""))
            elif entry.name.endswith(".tar"):
                tars.append(entry.path)

    def names(p: str) -> tuple[str, list[str]]:
        out = subprocess.run(["tar", "tf", p], capture_output=True, text=True, check=True)  # noqa: S603, S607
        return p, out.stdout.split()

    with ThreadPoolExecutor(8) as pool:
        for path, members in pool.map(names, tars):
            for mem in members:
                m = _MSA.search(mem)
                if m:
                    inv.setdefault(_norm(m.group(1)), ("tar", path, mem))
    return inv


def _cif_chain_seqs(path: Path) -> dict[str, str]:
    seq: dict[str, list[str]] = {}
    cols: list[str] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.split(".", 1)[1].strip())
            elif cols and line.startswith("ATOM"):
                d = dict(zip(cols, line.split(), strict=False))
                if d.get("label_atom_id") == "CA":
                    seq.setdefault(d["label_asym_id"], []).append(_AA.get(d["label_comp_id"], "X"))
    return {k: "".join(v) for k, v in seq.items()}


def plan() -> None:
    inv = _inventory()
    print(f"distributed MSA entities: {len(inv):,}", flush=True)

    roles: dict[str, set[str]] = {}
    homo = [line.split("\t", 1)[0] for line in (RAW / "homodimer" / "SOURCE.tsv").read_text().splitlines()[1:]]
    for e in homo:
        roles.setdefault(e, set()).add("homodimer")

    u2e: dict[str, list[str]] = {}
    with HOMO_META.open(newline="") as fh:
        for row in csv.DictReader(fh):
            u2e.setdefault(row["uniprotAccession"], []).append(_norm(row["modelEntityId"]))
    ambiguous = sum(1 for v in u2e.values() if len(v) > 1)
    print(f"UniProt accessions with >1 monomer entity: {ambiguous:,}", flush=True)
    by_seq = dict(line.split("\t") for line in SEQ_FOUND.read_text().splitlines())

    chain_rows = []
    hets = []
    with HET_META.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row["passes_quality_threshold"].strip().lower() != "true":
                continue
            het = row["modelEntityId"]
            hets.append(het)
            for chain, acc in (("A", row["uniprot_ac_1"]), ("B", row["uniprot_ac_2"])):
                cands = u2e.get(acc, [])
                if len(cands) == 1:
                    ent, how = cands[0], "uniprot"
                elif f"{het}:{chain}" in by_seq:
                    ent, how = by_seq[f"{het}:{chain}"], "sequence"
                elif cands:
                    ent, how = "|".join(cands), "uniprot_ambiguous"   # resolved by sequence in stage
                else:
                    ent, how = "", "unmapped"
                chain_rows.append((het, chain, acc, ent, how))
                for e in ent.split("|"):
                    if e:
                        roles.setdefault(e, set()).add("heterodimer")

    # chain sequences, for verifying every heterodimer chain's MSA by content
    het_dir = RAW / "heterodimer" / "cif"

    def seqs(het: str) -> list[tuple[str, str, str]]:
        n = het.split("-", 1)[1]
        return [(het, c, s) for c, s in _cif_chain_seqs(het_dir / n[-3:] / f"{het}-model_v1.cif.gz").items()]

    with ThreadPoolExecutor(16) as pool, CHAIN_SEQ.open("w") as out:
        for rows in pool.map(seqs, hets, chunksize=256):
            for r in rows:
                out.write("\t".join(r) + "\n")

    missing = [e for e in roles if e not in inv]
    with PLAN.open("w") as out:
        # sort by source so one array task owns each tar (a tar is then read once)
        for e in sorted(roles, key=lambda x: (inv.get(x, ("z", "", ""))[1], x)):
            kind, path, mem = inv.get(e, ("missing", "", ""))
            out.write(f"{e}\t{kind}\t{path}\t{mem}\t{','.join(sorted(roles[e]))}\n")
    with (WORK / "msa_stage_chainmap.tsv").open("w") as out:
        for r in chain_rows:
            out.write("\t".join(r) + "\n")
    print(f"planned {len(roles):,} MSAs  (not in the release: {len(missing):,} -> {missing[:5]})")
    print("chain links: " + ", ".join(f"{k} {sum(1 for r in chain_rows if r[4] == k):,}"
                                      for k in ("uniprot", "sequence", "uniprot_ambiguous", "unmapped")))


# ---------------------------------------------------------------- stage
def stage(task: int, ntasks: int) -> None:
    rows = [line.split("\t") for line in PLAN.read_text().splitlines()]
    # contiguous slices keep a tar's members together (the plan is sorted by source)
    lo, hi = task * len(rows) // ntasks, (task + 1) * len(rows) // ntasks
    mine = rows[lo:hi]
    linked = extracted = skipped = 0
    by_tar: dict[str, list[tuple[str, str]]] = {}
    for ent, kind, path, mem, _roles in mine:
        dest = _dest(ent)
        if dest.exists():
            skipped += 1
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        if kind == "dir":
            os.link(path, dest)
            linked += 1
        elif kind == "tar":
            by_tar.setdefault(path, []).append((ent, mem))
    scratch = Path(os.environ.get("TMPDIR", "/tmp")) / f"afm_msa_{task}"  # noqa: S108 (node-local)
    for tar, members in by_tar.items():
        scratch.mkdir(parents=True, exist_ok=True)
        lst = scratch / "members.txt"
        lst.write_text("".join(m + "\n" for _, m in members))
        subprocess.run(["tar", "xf", tar, "-C", str(scratch), "-T", str(lst)], check=True)  # noqa: S603, S607
        for ent, mem in members:
            shutil.move(str(scratch / mem), _dest(ent))
            extracted += 1
        shutil.rmtree(scratch, ignore_errors=True)
    print(f"task {task}: rows {len(mine):,}  linked {linked:,}  extracted {extracted:,}  already {skipped:,}")


# ---------------------------------------------------------------- finalize
def _query(path: Path) -> str:
    """Return the MSA's query sequence: the first record after any ``#`` header line.

    Reads only the head of the stream -- a deep MSA decompresses to tens of MB and only its
    first record matters here.
    """
    proc = subprocess.Popen(["zstd", "-dcq", str(path)], stdout=subprocess.PIPE)  # noqa: S603, S607
    head = proc.stdout.read(1 << 16).decode("ascii", "replace") if proc.stdout else ""
    proc.kill()
    proc.wait()
    lines = [ln for ln in head.split("\n") if ln and not ln.startswith("#")]
    return lines[1].strip() if len(lines) > 1 else ""


def finalize() -> None:
    # generated MSAs (entities the release does not ship): compressed alongside, in intermediate
    gen = INTER / "msa_generated"
    for a3m in sorted(gen.glob("*-msa_v1.a3m")):
        zst = a3m.with_name(a3m.name + ".zst")
        if not zst.exists():
            subprocess.run(["zstd", "-q", "-19", str(a3m), "-o", str(zst)], check=True)  # noqa: S603, S607
    generated = {p.name[: -len("-msa_v1.a3m.zst")]: p for p in gen.glob("*-msa_v1.a3m.zst")}

    def where(ent: str) -> Path:
        return generated.get(ent) or _dest(ent)

    plan_rows = [line.split("\t") for line in PLAN.read_text().splitlines()]
    have = sum(1 for r in plan_rows if _dest(r[0]).exists())  # release MSAs in raw
    print(f"staged {have:,} / planned {len(plan_rows):,}", flush=True)

    # verify every heterodimer chain link by content, resolving ambiguous UniProt links
    seq = {(h, c): s for h, c, s in (ln.split("\t") for ln in CHAIN_SEQ.read_text().splitlines())}
    chain_rows = [ln.split("\t") for ln in (WORK / "msa_stage_chainmap.tsv").read_text().splitlines()]

    def check(row: list[str]) -> list[str]:
        het, chain, acc, ent, how = row
        target = seq.get((het, chain), "")
        for cand in ent.split("|"):
            if cand and where(cand).exists() and _query(where(cand)) == target:
                return [het, chain, acc, cand, how if "|" not in ent else "uniprot+sequence", "ok"]
        return [het, chain, acc, ent, how, "MISMATCH"]

    with ThreadPoolExecutor(32) as pool:
        verified = list(pool.map(check, chain_rows, chunksize=512))
    bad = [r for r in verified if r[5] != "ok"]
    (INTER / "heterodimer").mkdir(parents=True, exist_ok=True)
    with (INTER / "heterodimer" / "chain_msa.tsv").open("w") as out:
        out.write("heterodimer_id\tchain\tuniprot\tmonomer_msa_entity\tlinked_by\tverified\n")
        for r in verified:
            out.write("\t".join(r) + "\n")
    print(f"heterodimer chains verified {len(verified) - len(bad):,} / {len(verified):,}  mismatch {len(bad):,}")

    with (RAW / "msa" / "SOURCE.tsv").open("w") as out:
        out.write("entity_id\tsource_kind\tsource_path\ttar_member\tused_by\n")
        for ent, kind, path, mem, roles in plan_rows:
            src_kind, src_path = kind, path
            if kind == "missing" and ent in generated:
                src_kind, src_path = "generated", str(generated[ent])
            out.write(f"{ent}\t{src_kind}\t{src_path}\t{mem}\t{roles}\n")

    n_gen = sum(1 for _ in gen.glob("*-msa_v1.a3m"))
    (RAW / "msa" / "README.md").write_text(
        f"""# msa — MSAs used by the AFM structures

{have:,} MSAs as the release ships them, one per monomer entity, in its own format
(`AF-<n>-msa_v1.a3m.zst`, one zstd frame per file), sharded on the last three digits:

    msa/<last 3 digits>/AF-<n>-msa_v1.a3m.zst
    SOURCE.tsv   entity, where it came from (release batch dir / tar member / generated), used by

A homodimer uses the MSA of its own entity. A heterodimer uses the MSAs of its two monomers
(the release concatenates them without pairing); the chain -> monomer link, every row
content-checked (MSA query == chain sequence), is
`BioMol/materials/intermediate/afdb_multimer/heterodimer/chain_msa.tsv`.

Source: the AFDB complex release msas/ (251124_1.2M, 251124_8.8M, 251209_13.4M) -- only the
entities used here, out of 23,425,744. Directory-stored files are hard links into the download
mirror; tar-only ones were extracted. {n_gen} entity(ies) the release does not ship have an MSA
we generated, in `BioMol/materials/intermediate/afdb_multimer/msa_generated/` (marked
`generated` in SOURCE.tsv).

The release's `pandemic_prep/` batch (200 tars, posted 2026-09-11) is not staged: no structure
here needs it.
""",
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    s = sub.add_parser("stage")
    s.add_argument("--task", type=int, required=True)
    s.add_argument("--ntasks", type=int, required=True)
    sub.add_parser("finalize")
    a = ap.parse_args()
    if a.cmd == "plan":
        plan()
    elif a.cmd == "stage":
        stage(a.task, a.ntasks)
    else:
        finalize()


if __name__ == "__main__":
    main()
