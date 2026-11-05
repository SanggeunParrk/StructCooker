"""Counts of the final pair table (asd/final/pairs_final.tsv): chains, clusters, antigen x heavy cluster combinations and their sizes.
Usage: asd_final_stats.py PAIRS_FINAL
"""

import collections
import sys

rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
ix = {k: i for i, k in enumerate(rows[0])}
rows = rows[1:]
col = lambda name: [r[ix[name]] for r in rows]
print(
    "pairs",
    len(rows),
    "| antigens",
    len(set(col("antigen_seq_id"))),
    "heavy",
    len(set(col("heavy_seq_id"))),
    "light",
    len({x for x in col("light_seq_id") if x}),
)
print(
    "with a light chain",
    sum(bool(x) for x in col("light_seq_id")),
    "| status",
    dict(collections.Counter(col("chain_status")).most_common()),
)
for name in ("cdr90", "cdr100"):
    h, l = col(f"heavy_cluster_{name}"), col(f"light_cluster_{name}")
    print(
        f"{name}: heavy clusters {len(set(h))}, light clusters {len({x for x in l if x})}, (heavy, light) {len(set(zip(h, l, strict=False)))}"
    )
a40, a100 = col("antigen_cluster40"), col("antigen_seq_id")
print("antigen clusters (40%)", len(set(a40)))
for name in ("cdr90", "cdr100"):
    h = col(f"heavy_cluster_{name}")
    for an, a in (("antigen 40%", a40), ("antigen 100%", a100)):
        m = collections.Counter(zip(a, h, strict=False))
        n = sorted(m.values())
        big = sum(v for v in n if v > 1000)
        print(
            f"{an} x heavy {name}: {len(m)} clusters; median {n[len(n) // 2]}, mean {sum(n) / len(n):.2f}, p99 {n[int(len(n) * 0.99)]}, max {n[-1]}; >1000: {sum(v > 1000 for v in n)} clusters hold {big} pairs; singletons {sum(v == 1 for v in n)}"
        )
m = collections.Counter(zip(a100, col("heavy_cluster_cdr100"), strict=False))
ds = collections.defaultdict(collections.Counter)
for r in rows:
    ds[(r[ix["antigen_seq_id"]], r[ix["heavy_cluster_cdr100"]])][r[ix["datasets"]]] += 1
print(
    "top antigen 100% x heavy cdr100:",
    [(c, dict(ds[k].most_common(2))) for k, c in m.most_common(6)],
)
