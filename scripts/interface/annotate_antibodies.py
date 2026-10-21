"""Classify every unique PDB protein sequence with ANARCI -> TSV cache (sbatch only).

usage: annotate_antibodies.py <out.tsv> [workers]
Input: $BIOMOL_ROOT/materials/intermediate/fasta/pdb_polypeptide_L.fasta. Output lines
``sequence<TAB>type`` (type empty = not an antibody / TCR); only lengths 70-600 are sent to ANARCI.
"""
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))


def classify(chunk):
    from anarci import anarci

    def one(batch):
        num, det, _ = anarci([(f"s{i}", s) for i, s in enumerate(batch)], scheme="imgt", output=False)
        return [(d[0].get("chain_type") if (n and d) else None) for n, d in zip(num, det, strict=True)]

    try:
        types = one(chunk)
    except Exception:  # noqa: BLE001  one bad sequence must not sink the chunk
        types = []
        for s in chunk:
            try:
                types += one([s])
            except Exception:  # noqa: BLE001
                types.append(None)
    return list(zip(chunk, types, strict=True))


def main():
    out, workers = Path(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 64
    seqs = set()
    with open(ROOT / "materials/intermediate/fasta/pdb_polypeptide_L.fasta") as f:
        for line in f:
            if not line.startswith(">"):
                s = line.strip()
                if 70 <= len(s) <= 600:
                    seqs.add(s)
    seqs = sorted(seqs)
    print(f"{len(seqs)} unique sequences", flush=True)
    chunks = [seqs[i : i + 200] for i in range(0, len(seqs), 200)]
    t0, n_ab = time.time(), 0
    with mp.get_context("fork").Pool(workers) as pool, open(out, "w") as fo:
        for i, res in enumerate(pool.imap_unordered(classify, chunks)):
            for s, t in res:
                fo.write(f"{s}\t{t or ''}\n")
                n_ab += bool(t)
            if i % 100 == 0:
                print(f"{i}/{len(chunks)} chunks, {n_ab} antibody/TCR, {time.time() - t0:.0f}s", flush=True)
    print(f"done: {n_ab} of {len(seqs)} classified as antibody/TCR, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
