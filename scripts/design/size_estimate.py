"""Stored size of the filter-passing items, from the source index (value_bytes) joined with items.tsv (sbatch only)."""

import sys

idx, items = sys.argv[1], sys.argv[2]
size = {}
with open(idx) as f:
    for ln in f:
        if ln[0] == "#" or ln.startswith("key\t"):
            continue
        k, b, _ = ln.rstrip("\n").split("\t")[:3]
        size[k] = int(b)
tot = ok = n = m = 0
by_len = {}
with open(items) as f:
    next(f)
    for ln in f:
        x = ln.rstrip("\n").split("\t")
        b = size[x[0]]
        tot += b
        if x[8] == "1":
            ok += b
            n += 1
            if x[7] == x[0]:
                m += b
                by_len[1] = by_len.get(1, 0) + 1
print(
    f"all {tot / 1e9:.0f} GB  passing {ok / 1e9:.0f} GB ({n} items, mean {ok / n / 1e3:.1f} KB)"
)
print(f"cluster representatives only (passing reps): {m / 1e9:.0f} GB")
