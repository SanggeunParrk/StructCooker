"""Mean local coordinates of a residue's own atoms in its backbone frame, for drawing (sbatch only).

usage: frames_skeleton.py <out.json> [n_entries=400]
Protein: N CA C O CB. Nucleic acid: sugar, phosphate and base atoms, purines (A, G) and pyrimidines (C, U, T) apart.
"""

import json
import os
import random
import sys
from pathlib import Path

import lmdb
import numpy as np

from structcooker.instructions.mols import CIFMolAttached
from structcooker.instructions.potential import frames as F
from structcooker.instructions.potential import pairs
from structcooker.instructions.transforms.codecs import from_bytes

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
DB = ROOT / "lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb"
SETS = {
    "protein": (
        ["N", "CA", "C", "O", "CB"],
        [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB")],
    ),
    "purine": (
        [
            "C1'",
            "C2'",
            "C3'",
            "C4'",
            "O4'",
            "C5'",
            "O5'",
            "O3'",
            "N9",
            "C4",
            "N3",
            "C2",
            "N1",
            "C6",
            "C5",
            "N7",
            "C8",
        ],
        [
            ("C1'", "C2'"),
            ("C2'", "C3'"),
            ("C3'", "C4'"),
            ("C4'", "O4'"),
            ("O4'", "C1'"),
            ("C4'", "C5'"),
            ("C5'", "O5'"),
            ("C3'", "O3'"),
            ("C1'", "N9"),
            ("N9", "C4"),
            ("C4", "N3"),
            ("N3", "C2"),
            ("C2", "N1"),
            ("N1", "C6"),
            ("C6", "C5"),
            ("C5", "C4"),
            ("C5", "N7"),
            ("N7", "C8"),
            ("C8", "N9"),
        ],
    ),
    "pyrimidine": (
        [
            "C1'",
            "C2'",
            "C3'",
            "C4'",
            "O4'",
            "C5'",
            "O5'",
            "O3'",
            "N1",
            "C2",
            "N3",
            "C4",
            "C5",
            "C6",
        ],
        [
            ("C1'", "C2'"),
            ("C2'", "C3'"),
            ("C3'", "C4'"),
            ("C4'", "O4'"),
            ("O4'", "C1'"),
            ("C4'", "C5'"),
            ("C5'", "O5'"),
            ("C3'", "O3'"),
            ("C1'", "N1"),
            ("N1", "C2"),
            ("C2", "N3"),
            ("N3", "C4"),
            ("C4", "C5"),
            ("C5", "C6"),
            ("C6", "N1"),
        ],
    ),
}
keys = [
    l.split("\t")[0]
    for l in open(str(DB) + ".index.tsv")
    if not l.startswith(("#", "key"))
]
random.Random(5).shuffle(keys)
env = lmdb.open(str(DB), readonly=True, lock=False)
acc = {k: {a: [] for a in v[0]} for k, v in SETS.items()}
for key in keys[: int(sys.argv[2]) if len(sys.argv) > 2 else 400]:
    with env.begin() as txn:
        rec = from_bytes(bytes(txn.get(key.encode())))
    ck = min(
        rec,
        key=lambda c: (
            int(c.split("_")[0]),
            int(c.split("_")[1]),
            c.split("_")[2] != ".",
            c,
        ),
    )
    cm = CIFMolAttached.from_dict(rec[ck]["cifmol_attached_dict"])
    rtype, poly = pairs.residue_types(cm)
    o, rot, ok = F.frames(cm, rtype)
    name, xyz = cm.atoms.id.value, cm.atoms.xyz.value
    ares = cm.index_table.atoms_to_residues(np.arange(len(name)))
    cls = np.full(len(rtype), "", dtype="U10")
    for r, t in enumerate(rtype.tolist()):
        if t < 0:
            continue
        label = pairs.TYPES[t]
        if t <= 20:
            cls[r] = "protein"
        elif label[1] in "AG":
            cls[r] = "purine"
        elif label[1] in "CUT":
            cls[r] = "pyrimidine"
    for k, (atoms, _) in SETS.items():
        for a in atoms:
            hit = np.flatnonzero(
                (name == a) & (cls[ares] == k) & ok[ares] & np.isfinite(xyz).all(axis=1)
            )
            if len(hit):
                acc[k][a].append(F.to_local(o[ares[hit]], rot[ares[hit]], xyz[hit]))
out = {}
for k, (atoms, bonds) in SETS.items():
    pos = {}
    for a in atoms:
        if acc[k][a]:
            v = np.concatenate(acc[k][a])
            pos[a] = [round(float(x), 3) for x in np.median(v, axis=0)] + [len(v)]
    out[k] = {"atoms": pos, "bonds": [b for b in bonds if b[0] in pos and b[1] in pos]}
json.dump(out, open(sys.argv[1], "w"))
for k in out:
    print(k, {a: v[:3] for a, v in list(out[k]["atoms"].items())[:6]})
