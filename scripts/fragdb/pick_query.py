"""Pick an antibody-antigen unit released in 2026 from the interface DB (sbatch only). usage: pick_query.py [index=0]"""
import csv
import hashlib
import os
import sys
from pathlib import Path

import numpy as np

from structcooker.interface.db import InterfaceDB

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")) / "lmdb/interface/pdb"
rows = []
for r in csv.DictReader(open(ROOT / "manifest_full.tsv"), delimiter="\t"):
    if not r["release_date"].startswith("2026") or r["has_ab"] != "1" or r["type"] != "3protein":
        continue
    ab = r["node_ab"].split("|")
    if "H" not in ab or not ({"K", "L"} & set(ab)) or ab.count("") != 1:
        continue
    if not 3000 <= int(r["n_atoms"]) <= 9000 or r["resolution"] in ("", "None", "nan") or float(r["resolution"]) > 3.5:
        continue
    rows.append(r)
rows.sort(key=lambda r: hashlib.md5(r["unit_id"].encode()).hexdigest())
print(f"{len(rows)} antibody-antigen units (H + L + antigen) released in 2026, resolution <= 3.5 A, 3,000-9,000 atoms")
db = InterfaceDB(ROOT)
for r in rows[:5]:
    rec = db.get(r["unit_id"], int(r["shard"]))
    nd = rec["nodes"]
    print(f"  {r['unit_id']}  released {r['release_date']}  {r['resolution']} A  atoms {r['n_atoms']}  interface residues {r['n_interface_res']}")
    print("      nodes:", [(nd["ab_class"][i] or "antigen", nd["chain_id"][i], nd["auth_asym_id"][i], len(nd["key"][i]) - 2) for i in range(len(nd["kind"]))])
