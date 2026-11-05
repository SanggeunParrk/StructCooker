"""Antibody sequence list ab_seq_id.tsv at 40% (sbatch only).

usage: make_ab_seq_id.py [--root $BIOMOL_ROOT]
An antibody sequence is a member of a ``cPDB_A...`` cluster of metadata/seq_cluster40_PDB.tsv (the release's antibody clusters: SabDab heavy / light
chains clustered on CDR3; see instructions/transforms/sequence.py antibody_cluster). Columns: seq_id, cluster40, anarci (H / K / L from
materials/intermediate/interface/antibody_cache.tsv, empty when ANARCI did not call it an antibody), length, sequence.
Written to materials/intermediate/interface/ab_seq_id.tsv.
"""

import argparse
import os
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"),
)
root = Path(ap.parse_args().root)

ab = {}
with open(root / "metadata/seq_cluster40_PDB.tsv") as f:
    for ln in f:
        c, m = ln.rstrip("\n").split("\t")
        if c.startswith("cPDB_A"):
            for s in m.split(","):
                ab[s] = c
print("antibody seq_ids", len(ab), "clusters", len(set(ab.values())))
anarci = {}
with open(root / "materials/intermediate/interface/antibody_cache.tsv") as f:
    for ln in f:
        x = ln.rstrip("\n").split("\t")
        if len(x) >= 2:
            anarci[x[0]] = x[1]
out = root / "materials/intermediate/interface/ab_seq_id.tsv"
tmp = out.with_suffix(".tmp")
rows = 0
with open(root / "metadata/seq_id_map.tsv") as f, open(tmp, "w") as o:
    o.write("seq_id\tcluster40\tanarci\tlength\tsequence\n")
    for ln in f:
        sid, seq = ln.rstrip("\n").split("\t")
        if sid in ab:
            o.write(f"{sid}\t{ab[sid]}\t{anarci.get(seq, '')}\t{len(seq)}\t{seq}\n")
            rows += 1
assert rows == len(ab), (rows, len(ab))
tmp.rename(out)
print("written", rows, out)
