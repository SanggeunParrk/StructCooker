"""MSA token alphabet as stored (schema E): letter -> token from the query row, gap token, over sampled records (sbatch only)."""

import sys
from collections import defaultdict

import lmdb
import numpy as np

from structcooker.instructions.transforms.codecs import from_bytes

for name, path, n in [("protein", sys.argv[1], 300), ("rna", sys.argv[2], 10_000)]:
    env = lmdb.open(path, readonly=True, lock=False, subdir=True)
    maps, allv = defaultdict(set), set()
    with env.begin() as t:
        for i, (_, v) in enumerate(t.cursor()):
            if i >= n:
                break
            s = from_bytes(bytes(v))["msa_dict"]["sequences"]
            q, a = (
                np.asarray(s["query_sequence"]).astype(str),
                np.asarray(s["aligned_sequences"]),
            )
            for letter, tok in zip(q.tolist(), a[0].tolist(), strict=False):
                maps[letter].add(tok)
            allv |= set(np.unique(a).tolist())
    print(name, "letter->token:", {k: sorted(v) for k, v in sorted(maps.items())})
    print(name, "token values seen:", sorted(allv))
