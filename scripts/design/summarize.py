"""Join per-item statistics with foldseek clusters and count what passes the filters (sbatch only).

usage: summarize.py DIR [--out DIR/items.tsv] [--all]
DIR holds part_*.tsv (build_stats) and fs/clusters.tsv (cluster_tar). Writes items.tsv: id, length, plddt, helix, strand, coil,
max_coil_run, cluster (representative id), pass (La-Proteina filters: pLDDT >= 80, 32 <= length <= 512, coil < 0.5, max coil run <= 20;
--all keeps every item, for sets that were filtered when they were made, e.g. Teddymer).
"""

import argparse
import glob
from collections import Counter
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("dir")
ap.add_argument("--out", default=None)
ap.add_argument("--all", action="store_true")
a = ap.parse_args()
d = Path(a.dir)
rep = {}
with open(d / "fs/clusters.tsv") as f:
    for ln in f:
        r, m = ln.rstrip("\n").split("\t")
        rep[m] = r
n = ok = miss = 0
reps, sizes = set(), Counter()
with open(a.out or d / "items.tsv", "w") as out:
    parts = sorted(glob.glob(str(d / "part_*.tsv")))
    with open(parts[0]) as f:
        has_asm = (
            len(f.readline().split("\t")) > 7
        )  # complexes by assembly: entry, asm follow the statistics
    out.write(
        "id\tlength\tplddt\thelix\tstrand\tcoil\tmax_coil_run\tcluster\tpass"
        + ("\tentry\tasm" if has_asm else "")
        + "\n"
    )
    for p in parts:
        with open(p) as f:
            for ln in f:
                x = ln.rstrip("\n").split("\t")
                c = rep.get(x[0])
                if c is None:
                    miss += 1
                    continue
                L, pl, co, mr = int(x[1]), float(x[2]), float(x[5]), int(x[6])
                keep = a.all or (pl >= 80 and 32 <= L <= 512 and co < 0.5 and mr <= 20)
                n += 1
                sizes[c] += 1
                if keep:
                    ok += 1
                    reps.add(c)
                out.write("\t".join(x[:7] + [c, str(int(keep))] + x[7:]) + "\n")
print(
    f"items {n}  no cluster {miss}  clusters {len(sizes)}  pass {ok}  clusters with a passing member {len(reps)}"
)
print(
    "cluster size: max",
    max(sizes.values()),
    " singletons",
    sum(v == 1 for v in sizes.values()),
)
if miss:
    raise SystemExit(f"{miss} items have no cluster; items.tsv is incomplete")
