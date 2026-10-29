"""RNA MSA token order in a FoldBench DB vs the release (sbatch only).

usage: rna_msa_token_order.py FOLDBENCH_RNA_LMDB RELEASE_RNA_LMDB SEQIDS.txt
For every key: letter -> token map from the query row (query_sequence vs aligned_sequences[0]); all records must agree,
and the FoldBench record must equal the release record byte for byte.
"""
import sys
from collections import Counter, defaultdict

import lmdb
import numpy as np

from structcooker.instructions.transforms.codecs import from_bytes, to_bytes

fb = lmdb.open(sys.argv[1], readonly=True, lock=False, subdir=True)
rel = lmdb.open(sys.argv[2], readonly=True, lock=False, subdir=True)
keys = [k.strip() for k in open(sys.argv[3]) if k.strip()]
maps = defaultdict(set)
orders = Counter()
c = Counter()
for k in keys:
    with fb.begin() as t:
        a = t.get(k.encode())
    with rel.begin() as t:
        b = t.get(k.encode())
    if a is None:
        c["missing_in_foldbench"] += 1
        continue
    c["found"] += 1
    c["bytes_equal_release"] += b is not None and bytes(a) == bytes(b)
    if b is not None and bytes(a) != bytes(b):
        c["decoded_equal_release"] += to_bytes(from_bytes(bytes(a))) == to_bytes(from_bytes(bytes(b)))
    s = from_bytes(bytes(a))["msa_dict"]["sequences"]
    q = np.asarray(s["query_sequence"]).astype(str)
    row0 = np.asarray(s["aligned_sequences"])[0]
    rec = {}
    for letter, tok in zip(q.tolist(), row0.tolist()):
        maps[letter].add(tok)
        rec[letter] = tok
    orders["".join(L for L, _ in sorted(((L, t) for L, t in rec.items() if L in "ACGU"), key=lambda x: x[1]))] += 1
print(dict(c))
print("letter -> token(s) over all records:", {k: sorted(v) for k, v in sorted(maps.items())})
print("records by letter order (token ascending):", dict(orders))
