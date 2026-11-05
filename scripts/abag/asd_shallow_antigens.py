"""Antigens with a shallow MSA: seq_id, depth, length, and what ASD calls them (target_name, uniprot, dataset, rows)."""

import collections
import glob
import sys

import pyarrow.parquet as pq

asd, depth_tsv, worklist, max_depth, out = (
    sys.argv[1],
    sys.argv[2],
    sys.argv[3],
    int(sys.argv[4]),
    sys.argv[5],
)
depth = {k: int(v) for k, v in (l.split() for l in open(depth_tsv))}
seq_of = dict(l.rstrip("\n").split("\t") for l in open(worklist))
want = {seq_of[k]: k for k, v in depth.items() if v <= max_depth}
info = {k: collections.defaultdict(collections.Counter) for k in want.values()}
for f in sorted(glob.glob(f"{asd}/part-*.parquet")):
    t = pq.read_table(f, columns=["dataset", "chains", "metadata"]).to_pylist()
    for r in t:
        s = "".join(
            (dict(r["chains"]["target_chains"] or []).get("full") or "").split()
        ).upper()
        k = want.get(s)
        if k is None:
            continue
        m = r["metadata"] or {}
        info[k]["name"][(m.get("target_name") or "").lower()] += 1
        info[k]["uniprot"][(m.get("target_uniprot") or "")] += 1
        info[k]["pdb"][(m.get("target_pdb") or "")] += 1
        info[k]["dataset"][r["dataset"]] += 1
with open(out, "w") as o:
    o.write(
        "seq_id\tdepth\tlength\trows\tdatasets\ttarget_name\tuniprot\tpdb\tsequence\n"
    )
    for s, k in sorted(want.items(), key=lambda x: (depth[x[1]], len(x[0]))):
        i = info[k]
        top = lambda c, n=2: "|".join(f"{a}" for a, _ in c.most_common(n))
        o.write(
            f"{k}\t{depth[k]}\t{len(s)}\t{sum(i['dataset'].values())}\t{top(i['dataset'], 4)}\t{top(i['name'])}\t{top(i['uniprot'])}\t{top(i['pdb'])}\t{s}\n"
        )
print(len(want))
