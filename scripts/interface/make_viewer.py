"""Build a self-contained 3D viewer page (3Dmol.js) from pick_examples.py output.

usage: make_viewer.py <examples.json> <out.html> [<db_root>]
The page embeds the example units (PDB text with the interface marker in the B-factor column) and
loads 3Dmol.js from cdnjs. Release numbers in the header are read from <db_root>/split_summary.md and
report.md when given.
"""
import json
import re
import sys
from pathlib import Path

examples, out = sys.argv[1], sys.argv[2]
root = Path(sys.argv[3]) if len(sys.argv) > 3 else None
units = json.load(open(examples))
head = {"units": "", "entries": "", "train": "", "valid": "", "test": ""}
if root and (root / "report.md").exists():
    m = re.search(r"Units: \*\*([\d,]+)\*\* from \*\*([\d,]+)\*\*", (root / "report.md").read_text())
    if m:
        head["units"], head["entries"] = m.group(1), m.group(2)
if root and (root / "split_summary.md").exists():
    for s in ("train", "valid", "test"):
        m = re.search(rf"\| {s} \| ([\d,]+) \|", (root / "split_summary.md").read_text())
        if m:
            head[s] = m.group(1)

TEMPLATE = (Path(__file__).parent / "viewer_template.html").read_text()
html = TEMPLATE.replace("/*DATA*/[]", json.dumps(units, separators=(",", ":"))).replace("/*HEAD*/{}", json.dumps(head))
Path(out).write_text(html)
print(f"{len(units)} units, {len(html) / 1e6:.1f} MB -> {out}")
