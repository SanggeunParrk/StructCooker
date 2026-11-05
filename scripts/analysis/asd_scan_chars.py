"""Count ASD chain sequences with a character outside A-Z, by chain, dataset and character."""

import collections
import glob
import re
import sys

import pyarrow.parquet as pq

root = sys.argv[1]
bad = re.compile(r"[^A-Z]")
n = collections.Counter()
chars = collections.Counter()
by_ds = collections.Counter()
ex = {}
for f in sorted(glob.glob(f"{root}/part-*.parquet")):
    t = pq.read_table(f, columns=["dataset", "chains"]).to_pylist()
    for r in t:
        for field, key in (
            ("antibody_chains", "heavy"),
            ("antibody_chains", "light"),
            ("target_chains", "full"),
        ):
            n[key, "rows"] += 1
            s = dict(r["chains"][field] or []).get(key)
            s = "".join((s or "").split())
            if not s:
                n[key, "empty"] += 1
                continue
            if s != s.upper():
                n[key, "lowercase"] += 1
            b = bad.findall(s.upper())
            if b:
                n[key, "bad"] += 1
                by_ds[key, r["dataset"]] += 1
                chars.update(set(b))
                ex.setdefault((key, b[0]), s[:80])
for k, v in sorted(n.items()):
    print(k, v)
print("chars", chars.most_common())
print("by dataset", by_ds.most_common(20))
for k, v in ex.items():
    print(k, v)
