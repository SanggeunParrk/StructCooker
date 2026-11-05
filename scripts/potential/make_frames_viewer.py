"""Build the self-contained frame-contact viewer page from frames_reduce.py output (one file per contact definition).

usage: make_frames_viewer.py <sidechain.json> <all.json> <out.html>
"""

import json
import sys
from pathlib import Path

data = json.dumps(
    {"sc": json.load(open(sys.argv[1])), "all": json.load(open(sys.argv[2]))},
    separators=(",", ":"),
)
html = (
    (Path(__file__).parent / "frames_viewer_template.html")
    .read_text()
    .replace("/*DATA*/null", data)
)
Path(sys.argv[3]).write_text(html)
print(f"{len(html) / 1e6:.1f} MB -> {sys.argv[3]}")
