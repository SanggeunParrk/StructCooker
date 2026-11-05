"""Clusters of the ASD pair table (asd/pairs/pairs_kd1um.tsv, recipe asd/pairs_kd1um): how many antigen / heavy / light
clusters, their combinations, and the members per (antigen cluster, heavy cluster). Reads only the pair table and the two
cluster files, so a new clustering is one run of this script.

Usage: asd_pair_cluster_stats.py PAIRS CLUSTER_ASD|seq_id CLUSTER_ASB [MEMBERS_OUT]
``seq_id`` as the antigen cluster file: an antigen is its own cluster (100% identity).
A few antibody ids sit in two clusters (chains ANARCI could not number); the first cluster is kept.
"""

import collections
import sys

pairs, clu_asd, clu_asb = sys.argv[1:4]
members_out = sys.argv[4] if len(sys.argv) > 4 else None


def clusters(path):
    m = {}
    for line in open(path):
        c, members = line.rstrip("\n").split("\t")
        for x in members.split(","):
            m.setdefault(x, c)
    return m


ca = {} if clu_asd == "seq_id" else clusters(clu_asd)
cb = clusters(clu_asb)
antigen_of = (lambda a: a) if clu_asd == "seq_id" else (lambda a: ca[a])
rows = []
datasets = collections.defaultdict(collections.Counter)
for line in list(open(pairs))[1:]:
    a, h, l, ds = line.rstrip("\n").split("\t")
    r = (antigen_of(a), cb[h], cb[l] if l else "")
    rows.append(r)
    datasets[r[:2]][ds] += 1
print("pairs", len(rows), " distinct triples", len(set(rows)))
print(
    "antigen clusters:",
    len({r[0] for r in rows}),
    " heavy clusters:",
    len({r[1] for r in rows}),
    " light clusters:",
    len({r[2] for r in rows if r[2]}),
)
print(
    "antibody (heavy, light) cluster combinations:", len({(r[1], r[2]) for r in rows})
)
m = collections.Counter((r[0], r[1]) for r in rows)
print("antigen cluster x heavy cluster:", len(m), " x light cluster:", len(set(rows)))
if members_out:
    with open(members_out, "w") as o:
        for (a, h), n in sorted(m.items()):
            o.write(f"{a}\t{h}\t{n}\n")
n = sorted(m.values())
print(
    "members per antigen x heavy cluster: median",
    n[len(n) // 2],
    "mean",
    round(sum(n) / len(n), 2),
    "p90",
    n[int(len(n) * 0.9)],
    "p99",
    n[int(len(n) * 0.99)],
    "max",
    n[-1],
)
edges = [1, 2, 3, 4, 6, 11, 21, 51, 101, 501, 1001, 10**9]
labels = [
    "1",
    "2",
    "3",
    "4-5",
    "6-10",
    "11-20",
    "21-50",
    "51-100",
    "101-500",
    "501-1000",
    ">1000",
]
for i, lb in enumerate(labels):
    c = [x for x in n if edges[i] <= x < edges[i + 1]]
    print(f"  {lb:>8} clusters={len(c):6d} pairs={sum(c):7d}")
print("top clusters:")
for (a, h), c in m.most_common(6):
    print("  ", c, a, h, dict(datasets[(a, h)].most_common(3)))
