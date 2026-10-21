"""Pick representative units and convert them to PDB text for a 3D viewer (sbatch only).

usage: pick_examples.py <db_root> <out.json>
PDB text conventions: one chain letter per node (waters: W), residue number = position in the node's
sequence, B-factor column = marker (100 interface residue, 50 other contact atom, 0 rest; water 20, bridging 80).
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

from structcooker.interface.db import InterfaceDB

root, out = Path(sys.argv[1]), sys.argv[2]
rows = list(csv.DictReader(open(root / "manifest_full.tsv"), delimiter="\t"))
split = {}
if (root / "split.tsv").exists():
    for r in csv.DictReader(open(root / "split.tsv"), delimiter="\t"):
        split[r["unit_id"]] = (r["split"], r["novel_protein"], r["novel_pair"])
by_id = {r["unit_id"]: r for r in rows}
f = lambda r, k: float(r[k]) if r[k] not in ("", "None") else 99.0  # noqa: E731
base = [r for r in rows if 1200 <= int(r["n_atoms"]) <= 9000 and r["truncated"] == "False" and f(r, "resolution") <= 2.6]


def h(r):
    return hashlib.md5(r["unit_id"].encode()).hexdigest()  # noqa: S324


CATS = [
    ("antibody - antigen", lambda r: r["has_ab"] == "1" and r["n_protein"] == "3" and r["n_ligand"] == "0" and r["n_na"] == "0" and "H" in r["node_ab"].split("|")),
    ("TCR - pMHC", lambda r: {"A", "B"} <= set(r["node_ab"].split("|")) and int(r["n_protein"]) >= 4),
    ("protein - dsDNA + metal ion", lambda r: r["n_na"] == "2" and int(r["n_duplex_pairs"]) >= 1 and int(r["n_ion"]) >= 1 and int(r["n_protein"]) <= 2),
    ("protein - dsDNA", lambda r: r["n_na"] == "2" and int(r["n_duplex_pairs"]) >= 1 and r["n_ion"] == "0" and int(r["n_protein"]) <= 2),
    ("protein - single-strand NA", lambda r: r["n_na"] == "1" and int(r["n_protein"]) <= 2 and r["n_ligand"] == "0"),
    ("protein + ATP / ADP", lambda r: any(c in ("ATP", "ADP", "ANP", "AGS") for c in r["node_ccd"].split("|")) and int(r["n_protein"]) <= 2),
    ("protein + heme", lambda r: "HEM" in r["node_ccd"].split("|") and int(r["n_protein"]) <= 3),
    ("metal ion bridging two chains", lambda r: r["n_protein"] == "2" and int(r["n_ion"]) >= 1 and r["n_ligand"] == "0" and r["n_na"] == "0"),
    ("glycan on protein", lambda r: any("-" in c for c in r["node_ccd"].split("|")) and int(r["n_protein"]) <= 2),
    ("homodimer with bridging waters", lambda r: r["type"] == "2protein" and int(r["n_bridging_water"]) >= 15 and f(r, "resolution") <= 1.8),
    ("homotrimer", lambda r: r["type"] == "3protein" and int(r["n_water"]) >= 100),
]
MUST = ["3hfm_1_1_0", "1tsr_1_1_0"]
picked = []
for name, pred in CATS:
    cand = sorted((r for r in base if pred(r) and all(r["pdb_id"] != p["pdb_id"] for _, p in picked)), key=h)
    if cand:
        picked.append((name, cand[0]))
for uid in MUST:
    if uid in by_id and all(uid != r["unit_id"] for _, r in picked):
        picked.append(("reference structure", by_id[uid]))

LETTERS = "ABCDEFGHIJKLMNOPQRSTUVXYZ"
db = InterfaceDB(root)
units = []
for name, r in picked:
    rec = db.get(r["unit_id"], int(r["shard"]))
    meta, nodes, at, rs, w = rec["meta"], rec["nodes"], rec["atoms"], rec["residues"], rec["water"]
    lines, n = [], 0
    for i in range(len(at["xyz"])):
        ri, nd = int(at["res"][i]), int(at["node"][i])
        kind = nodes["kind"][nd]
        name_a, el = str(at["atom_name"][i]), str(at["element"][i])
        an = f" {name_a:<3}" if len(el) == 1 and len(name_a) < 4 else f"{name_a:<4}"
        b = 100.0 if rs["interface"][ri] else (50.0 if at["contact"][i] else 0.0)
        x, y, z = at["xyz"][i]
        n += 1
        lines.append(f"{'ATOM  ' if kind in ('protein', 'na') else 'HETATM'}{n % 100000:5d} {an} {rs['chem_comp_id'][ri]!s:>3} "
                     f"{LETTERS[nd % 25]}{int(rs['cif_idx'][ri]) % 10000:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{b:6.2f}          {el:>2}")
    for j in range(len(w["xyz"])):
        x, y, z = w["xyz"][j]
        n += 1
        lines.append(f"HETATM{n % 100000:5d}  O   HOH W{(j + 1) % 10000:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{80.0 if w['bridging'][j] else 20.0:6.2f}           O")
    sp = split.get(r["unit_id"], ("", "", ""))
    units.append({
        "id": r["unit_id"], "category": name, "type": r["type"], "split": sp[0], "novel_protein": sp[1], "novel_pair": sp[2],
        "resolution": meta["resolution"], "release_date": meta["release_date"],
        "n_atoms": len(at["xyz"]), "n_interface_res": int(r["n_interface_res"]), "n_water": len(w["xyz"]),
        "n_bridging": int(w["bridging"].sum()), "n_duplex": int(r["n_duplex_pairs"]),
        "nodes": [{"chain": LETTERS[i % 25], "kind": nodes["kind"][i], "auth": nodes["auth_asym_id"][i],
                   "label": (nodes["key"][i][2:] if nodes["kind"][i] in ("ligand", "ion") else f"{len(nodes['key'][i]) - 2} residues"),
                   "ab": nodes["ab_class"][i]} for i in range(len(nodes["kind"]))],
        "pdb": "\n".join(lines) + "\nEND\n",
    })
    print(r["unit_id"], name, r["type"], len(at["xyz"]), "atoms")
json.dump(units, open(out, "w"))
print(len(units), "units ->", out)
