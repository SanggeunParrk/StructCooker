"""Depth (rows of each a3m, query included) of every <seq_id>.a3m under a tree -> TSV seq_id, depth."""

import sys
from pathlib import Path

root, out = Path(sys.argv[1]), Path(sys.argv[2])
rows = []
for p in sorted(root.rglob("*.a3m")):
    with p.open("rb") as f:
        rows.append((p.stem, sum(1 for line in f if line.startswith(b">"))))
out.write_text("".join(f"{k}\t{v}\n" for k, v in rows))
print(len(rows))
