"""S1 -> S2 -> S3 funnel of ASD v2.0 (the SelfDistillation definitions), with and without the antigen filter.

S1: metadata.target_pdb empty (no experimental structure). S2: a heavy chain and an antigen sequence are present
(and valid, A-Z). S3: distinct (heavy, light, antigen). The antigen filter: the antigen is on the work list and is not
(MSA depth <= max_depth and no template). A final 'light only' filter is applied last to show it removes nothing S2 left.
"""
import collections, glob, re, sys

import pyarrow.parquet as pq

asd, worklist, depth_tsv, tmpl_tsv, max_depth = sys.argv[1:6]
max_depth = int(max_depth)
seq_of = dict(l.rstrip("\n").split("\t") for l in open(worklist))
depth = {k: int(v) for k, v in (l.split() for l in open(depth_tsv))}
tmpl = {k: int(v) for k, v in (l.split() for l in open(tmpl_tsv))}
keep = {k for k in seq_of if not (depth[k] <= max_depth and tmpl.get(k, 0) == 0)}
id_of = {s: k for k, s in seq_of.items()}
ok = re.compile(r"[A-Z]+")
norm = lambda x: "".join((x or "").split()).upper()
n = collections.Counter(); t_all = {}; t_keep = {}; ds_keep = collections.Counter()
for f in sorted(glob.glob(f"{asd}/part-*.parquet")):
    for r in pq.read_table(f, columns=["dataset", "chains", "metadata"]).to_pylist():
        n["S0"] += 1
        if (r["metadata"] or {}).get("target_pdb"):
            continue
        n["S1"] += 1
        ab = dict(r["chains"]["antibody_chains"] or []); tg = dict(r["chains"]["target_chains"] or [])
        h, l, a = norm(ab.get("heavy")), norm(ab.get("light")), norm(tg.get("full"))
        if not h or not a or not ok.fullmatch(h) or not ok.fullmatch(a) or (l and not ok.fullmatch(l)):
            continue
        n["S2"] += 1
        t_all.setdefault((h, l, a), r["dataset"])
        k = id_of.get(a)
        if k is None or k not in keep:
            continue
        n["S2+antigen"] += 1
        t_keep.setdefault((h, l, a), r["dataset"])
print("rows   S0", n["S0"], " S1", n["S1"], " S2", n["S2"], " S2 + antigen filter", n["S2+antigen"])
print("S3 distinct pairs: all antigens", len(t_all), "| after the antigen filter", len(t_keep))
lo = sum(1 for (h, l, a) in t_keep if not h)
print("light-only among them (filter removes):", lo, "-> final", len(t_keep) - lo)
kinds = collections.Counter("heavy+light" if l else "heavy only" for (h, l, a) in t_keep)
print("final composition:", dict(kinds))
print("final by dataset:", dict(collections.Counter(t_keep.values()).most_common()))
print("final distinct heavy", len({h for h, l, a in t_keep}), "light", len({l for h, l, a in t_keep if l}), "antigen", len({a for h, l, a in t_keep}))
