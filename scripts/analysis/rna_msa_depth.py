"""Depth of the PDB RNA alignments: ours (materials/reference/msa/R, nhmmer a3m) vs OpenFold3's
(s3://openfold3-data/pdb_training_set/alignment_arrays, nucleotide_collection / rfam / rnacentral hits).

Mapping: an OpenFold file is ``<pdb>_<chain>.npz``. The chain is looked up in ``cif_pdb.fasta`` by its
label asym id (then its author id), keeping candidates of type polyribonucleotide whose length equals
the width of the OpenFold ``msa`` array; the query sequence, sanitised the way our nhmmer input was
(letters other than A C G U N X -> X), gives our ``R`` seq_id through the a3m queries.

Depth = number of sequences in the alignment. Ours counts the a3m records (query included). OpenFold's is
the rows of the hit arrays, summed over sources ("rows") and as distinct sequences ("unique").

The OpenFold npz are pickled object arrays from the official bucket; the build reads them with
``allow_pickle=True`` too.

    python rna_msa_depth.py --out <dir> [--workers 16]
"""
from __future__ import annotations

import argparse
import glob
import os
import zipfile
from multiprocessing import get_context
from pathlib import Path

import numpy as np

ROOT = Path("/data/shared/cssb_data/BioMol_clean/materials")
OPENFOLD = ROOT / "raw/rna_alignment_arrays"
OURS = ROOT / "reference/msa/R"
FASTA = ROOT / "intermediate/fasta/cif_pdb.fasta"
RNA_SOURCES = ("nucleotide_collection_hits", "rfam_hits", "rnacentral_hits")
ALLOWED = set("ACGUNX")


def sanitize(seq: str) -> str:
    return "".join(c if c in ALLOWED else "X" for c in seq)


def load_ours() -> tuple[dict[str, str], dict[str, dict[str, int]]]:
    """sanitised query -> R seq_id, and per R id the a3m depth and its sources."""
    by_seq: dict[str, str] = {}
    depth: dict[str, dict[str, int]] = {}
    for path in glob.glob(str(OURS / "*/*/*.a3m")):
        rid = os.path.basename(path)[:-4]
        n = rnac = 0
        with open(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    n += 1
                    rnac += line.startswith(">URS")
        fa = path[:-4] + ".fasta"
        if os.path.exists(fa):
            by_seq[open(fa).read().split("\n")[1].strip()] = rid
        depth[rid] = {"depth": n, "rnacentral": rnac, "other": n - rnac}
    return by_seq, depth


def load_fasta() -> tuple[dict, dict]:
    label: dict[tuple[str, str], list[tuple[str, str]]] = {}
    auth: dict[tuple[str, str], list[tuple[str, str]]] = {}
    with open(FASTA) as handle:
        hdr = None
        for line in handle:
            if line.startswith(">"):
                hdr = line[1:].strip()
            elif hdr:
                idp, typ, a = (x.strip() for x in hdr.split("|"))
                pdb, chain, _alt = idp.split("_", 2)
                if typ == "polyribonucleotide":
                    label.setdefault((pdb.lower(), chain), []).append(line.strip())
                    auth.setdefault((pdb.lower(), a.removeprefix("Auth:")), []).append(line.strip())
                hdr = None
    return label, auth


_G: dict = {}


def _init(label, auth, by_seq):
    _G.update(label=label, auth=auth, by_seq=by_seq)


def one(path: str) -> dict | None:
    stem = os.path.basename(path)[:-4]
    with zipfile.ZipFile(path) as z:
        names = {i.filename.removesuffix(".npy") for i in z.infolist()}
    if not names & set(RNA_SOURCES):
        return None
    pdb, chain = stem.split("_", 1)
    z = np.load(path, allow_pickle=True)
    counts: dict[str, int] = {}
    uniq: set[str] = set()
    width = None
    for k in RNA_SOURCES:
        if k not in z.files:
            counts[k] = 0
            continue
        d = z[k].item()
        m = d["msa"]
        counts[k] = int(m.shape[0])
        if m.shape[0]:
            width = int(m.shape[1])
            uniq.update("".join(r) for r in m)
    cands = _G["label"].get((pdb, chain), []) or _G["auth"].get((pdb, chain), [])
    seqs = [s for s in cands if width is not None and len(s) == width]
    rid = _G["by_seq"].get(sanitize(seqs[0])) if seqs else None
    return {"file": stem, "width": width, "rid": rid, "mapped_rna": bool(cands), **{f"n_{k}": v for k, v in counts.items()},
            "rows": sum(counts.values()), "unique": len(uniq)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0, help="only the first N files (smoke test)")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    by_seq, ours = load_ours()
    label, auth = load_fasta()
    print(f"ours: {len(ours):,} RNA alignments; {len(by_seq):,} distinct queries", flush=True)
    files = sorted(glob.glob(str(OPENFOLD / "*.npz")))
    if args.limit:
        files = files[: args.limit]
    rows = []
    with get_context("forkserver").Pool(args.workers, initializer=_init, initargs=(label, auth, by_seq)) as pool:
        for i, r in enumerate(pool.imap_unordered(one, files, chunksize=8), 1):
            if r:
                rows.append(r)
            if i % 2000 == 0:
                print(f"  {i}/{len(files)}", flush=True)
    import csv
    with open(args.out / "per_chain.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["openfold_file", "width", "seq_id", "of_nt", "of_rfam", "of_rnacentral", "of_rows", "of_unique",
                    "our_depth", "our_rnacentral", "our_other"])
        for r in rows:
            o = ours.get(r["rid"], {})
            w.writerow([r["file"], r["width"], r["rid"] or "", r["n_nucleotide_collection_hits"], r["n_rfam_hits"],
                        r["n_rnacentral_hits"], r["rows"], r["unique"], o.get("depth", ""), o.get("rnacentral", ""), o.get("other", "")])
    print(f"done: {len(rows):,} OpenFold RNA files; mapped to ours: {sum(1 for r in rows if r['rid'] in ours):,}")


if __name__ == "__main__":
    main()
