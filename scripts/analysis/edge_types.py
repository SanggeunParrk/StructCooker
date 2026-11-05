"""What each MiniWorld sampler type is made of: cluster-letter combinations per edge type (sbatch only).

usage: edge_types.py train_item.tsv
"""

import re
import sys
from collections import Counter

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from expected_draws import edge_type

edges = Counter()
with open(sys.argv[1]) as f:
    next(f)
    for ln in f:
        x = ln.split("\t", 3)
        edges[x[1] if x[2] == "None" else f"{x[1]}_{x[2]}"] += 1
rows, clusters = Counter(), Counter()
for e, n in edges.items():
    combo = "+".join(sorted(re.findall(r"c(?:PDB_)?([A-Z])\d", e)))
    key = (edge_type(e), combo)
    rows[key] += n
    clusters[key] += 1
print("type\tletters\trows\tclusters")
for (t, combo), n in sorted(rows.items(), key=lambda kv: (kv[0][0], -kv[1])):
    print(f"{t}\t{combo}\t{n:,}\t{clusters[(t, combo)]:,}")
