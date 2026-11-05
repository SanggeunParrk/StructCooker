"""What ASD calls the antigens of a seq_id list: seq_id, length, rows, datasets, target_name, uniprot, pdb, sequence."""
import collections, glob, sys

import pyarrow.parquet as pq

asd, worklist, idlist, out = sys.argv[1:5]
seq_of = dict(l.rstrip("\n").split("\t") for l in open(worklist))
ids = {l.strip() for l in open(idlist) if l.strip()}
want = {seq_of[k]: k for k in ids}
info = {k: collections.defaultdict(collections.Counter) for k in ids}
for f in sorted(glob.glob(f"{asd}/part-*.parquet")):
    for r in pq.read_table(f, columns=["dataset", "chains", "metadata"]).to_pylist():
        k = want.get("".join((dict(r["chains"]["target_chains"] or []).get("full") or "").split()).upper())
        if k is None:
            continue
        m = r["metadata"] or {}
        info[k]["name"][(m.get("target_name") or "").lower()] += 1
        info[k]["uniprot"][m.get("target_uniprot") or ""] += 1
        info[k]["pdb"][m.get("target_pdb") or ""] += 1
        info[k]["dataset"][r["dataset"]] += 1
top = lambda c, n=2: "|".join(a for a, _ in c.most_common(n))
with open(out, "w") as o:
    o.write("seq_id\tlength\trows\tdatasets\ttarget_name\tuniprot\tpdb\tsequence\n")
    for s, k in sorted(want.items(), key=lambda x: (-sum(info[x[1]]["dataset"].values()), len(x[0]))):
        i = info[k]
        o.write(f"{k}\t{len(s)}\t{sum(i['dataset'].values())}\t{top(i['dataset'], 4)}\t{top(i['name'])}\t{top(i['uniprot'])}\t{top(i['pdb'])}\t{s}\n")
