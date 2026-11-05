"""Keep the antigens that are NOT (MSA depth <= max_depth AND no template); count the antigens and the antigen-antibody pairs.

A pair is a distinct (heavy, light, antigen) triple of ASD rows whose antigen is kept. A chain that is empty is
absent (a nanobody has no light chain); a row with a chain ASD ships with a character outside A-Z is dropped, since
it has no seq_id (asd/antibody_fasta skips it).
"""

import collections
import glob
import re
import sys

import pyarrow.parquet as pq

asd, worklist, depth_tsv, tmpl_tsv, max_depth, out_keep = sys.argv[1:7]
max_depth = int(max_depth)
seq_of = dict(l.rstrip("\n").split("\t") for l in open(worklist))
depth = {k: int(v) for k, v in (l.split() for l in open(depth_tsv))}
tmpl = {k: int(v) for k, v in (l.split() for l in open(tmpl_tsv))}
drop = {k for k in seq_of if depth[k] <= max_depth and tmpl.get(k, 0) == 0}
keep = {k for k in seq_of if k not in drop}
print(
    f"worklist antigens {len(seq_of)}  dropped (depth<={max_depth} and 0 templates) {len(drop)}  kept {len(keep)}"
)
print(
    f"  of the dropped: depth<={max_depth} {sum(depth[k] <= max_depth for k in seq_of)}, 0 templates {sum(tmpl.get(k, 0) == 0 for k in seq_of)}"
)
open(out_keep, "w").write("".join(f"{k}\n" for k in sorted(keep)))
id_of = {s: k for k, s in seq_of.items()}
letters = re.compile(r"[A-Z]+")
norm = lambda x: "".join((x or "").split()).upper()
rows = collections.Counter()
triples = {}
bad = 0
outside = collections.Counter()
for f in sorted(glob.glob(f"{asd}/part-*.parquet")):
    for r in pq.read_table(f, columns=["dataset", "chains"]).to_pylist():
        ab = dict(r["chains"]["antibody_chains"] or [])
        tg = dict(r["chains"]["target_chains"] or [])
        h, l, a = norm(ab.get("heavy")), norm(ab.get("light")), norm(tg.get("full"))
        rows["all"] += 1
        k = id_of.get(a)
        if k is None:
            outside["antigen not on the worklist (<10 or >3000 aa, invalid)"] += 1
            continue
        if k not in keep:
            outside["antigen dropped (shallow MSA, no template)"] += 1
            continue
        if (h and not letters.fullmatch(h)) or (l and not letters.fullmatch(l)):
            bad += 1
            continue
        kind = (
            "heavy+light"
            if h and l
            else "heavy only"
            if h
            else "light only"
            if l
            else "no antibody"
        )
        rows[kind] += 1
        triples.setdefault((h, l, k), (kind, r["dataset"]))
print("\nASD rows", rows["all"], "->", dict(outside), "invalid antibody chain", bad)
print(
    "rows with a kept antigen:",
    {k: v for k, v in rows.items() if k != "all"},
    " total",
    sum(v for k, v in rows.items() if k != "all"),
)
kinds = collections.Counter(v[0] for v in triples.values())
print("\ndistinct (heavy, light, antigen) pairs:", len(triples))
for k, v in kinds.most_common():
    print(f"  {k:12s} {v:8d}")
print(
    "  distinct heavy chains",
    len({t[0] for t in triples if t[0]}),
    " distinct light",
    len({t[1] for t in triples if t[1]}),
    " distinct antigens",
    len({t[2] for t in triples}),
)
ds = collections.Counter(v[1] for v in triples.values())
print("by dataset:", ", ".join(f"{k} {v}" for k, v in ds.most_common(10)))
