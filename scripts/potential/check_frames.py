"""Check the backbone frames: a residue's own CB (protein) / C2 or C4 (nucleic acid) must sit at the same local
coordinates in every frame (sbatch only). usage: check_frames.py [n_entries=60]"""
import os
import random
import sys
from pathlib import Path

import lmdb
import numpy as np

from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.mols import CIFMolAttached
from structcooker.potential import frames as F
from structcooker.potential import pairs

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
DB = ROOT / "lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb"
keys = [l.split("\t")[0] for l in open(str(DB) + ".index.tsv") if not l.startswith(("#", "key"))]
random.Random(3).shuffle(keys)
env = lmdb.open(str(DB), readonly=True, lock=False)
acc = {t: [] for t in range(pairs.T)}
for key in keys[: int(sys.argv[1]) if len(sys.argv) > 1 else 60]:
    with env.begin() as txn:
        rec = from_bytes(bytes(txn.get(key.encode())))
    ck = min(rec, key=lambda c: (int(c.split("_")[0]), int(c.split("_")[1]), c.split("_")[2] != ".", c))
    cm = CIFMolAttached.from_dict(rec[ck]["cifmol_attached_dict"])
    rtype, poly = pairs.residue_types(cm)
    rep = pairs.representative_atoms(cm, rtype)
    o, rot, ok = F.frames(cm, rtype)
    xyz = cm.atoms.xyz.value
    good = np.flatnonzero(ok & (rep >= 0))
    loc = F.to_local(o[good], rot[good], xyz[rep[good]])
    # proper rotations only: det = +1
    assert np.allclose(np.linalg.det(rot[good]), 1.0, atol=1e-5), "a frame is a reflection"
    for t, p in zip(rtype[good], loc):
        acc[int(t)].append(p)
print("own representative atom in the frame of its residue (x, y, z in A): mean +- sd over residues")
for t in (0, 2, 5, 9, 19, 20, 21, 22, 23, 24, 26, 27, 28, 29):
    a = np.array(acc[t])
    if len(a) > 20:
        print(f"  {pairs.TYPES[t]:>3s} n={len(a):6d}  " + "  ".join(f"{m:6.2f}+-{s:.2f}" for m, s in zip(a.mean(0), a.std(0))))
