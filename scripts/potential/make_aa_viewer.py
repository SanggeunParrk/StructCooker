"""Build the self-contained amino-acid pair distribution page from analyze_hist.py output of both definitions.

usage: make_aa_viewer.py <A.json (5 A contact)> <B.json (all pairs within 15 A)> <out.html>
"""
import json
import sys
from pathlib import Path

data = json.dumps({"A": json.load(open(sys.argv[1])), "B": json.load(open(sys.argv[2]))}, separators=(",", ":"))
html = (Path(__file__).parent / "aa_viewer_template.html").read_text().replace("/*DATA*/null", data)
Path(sys.argv[3]).write_text(html)
print(f"{len(html) / 1e3:.0f} KB -> {sys.argv[3]}")
