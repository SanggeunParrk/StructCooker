"""Build the protein-NA / NA-NA page from analyze_na.py output of both definitions.

usage: make_na_viewer.py <A.json> <B.json> <out.html>
"""
import json
import sys
from pathlib import Path

data = json.dumps({"A": json.load(open(sys.argv[1])), "B": json.load(open(sys.argv[2]))}, separators=(",", ":"))
html = (Path(__file__).parent / "na_viewer_template.html").read_text().replace("/*DATA*/null", data)
Path(sys.argv[3]).write_text(html)
print(f"{len(html) / 1e3:.0f} KB -> {sys.argv[3]}")
