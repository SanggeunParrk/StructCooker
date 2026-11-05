"""Check a built design DB against its source and against a fresh secondary-structure computation (sbatch only).

usage: validate_db.py DIR NAME SRC [--n 500]
"""
import argparse
import json
import random
import sys
from pathlib import Path

import lmdb
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from build_stats import shard_paths  # noqa: E402
from structcooker.instructions.design.ss import secondary_structure  # noqa: E402
from structcooker.instructions.transforms.codecs import from_bytes  # noqa: E402
from structcooker.instructions.mols import CIFMolAttached  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("dir")
ap.add_argument("name")
ap.add_argument("src")
ap.add_argument("--n", type=int, default=500)
a = ap.parse_args()
coll = Path(a.dir) / f"{a.name}.lmdb"
meta = json.loads((coll / "shards.json").read_text())
envs = [lmdb.open(str(coll / s["path"]), readonly=True, lock=False, subdir=True) for s in meta["shards"]]
src = [lmdb.open(p, readonly=True, lock=False, subdir=Path(p).is_dir()) for p in shard_paths(a.src)]
keys = []
for e in envs:
    with e.begin() as t:
        keys += [(k, e) for k in t.cursor().iternext(values=False)]
assert len(keys) == meta["entries"], (len(keys), meta["entries"])
manifest = {}
for ln in open(Path(a.dir) / "manifest.tsv").read().splitlines()[1:]:
    x = ln.split("\t")
    manifest[x[0]] = x
assert len(manifest) == len(keys), (len(manifest), len(keys))
random.Random(0).shuffle(keys)
bad = 0
for k, e in keys[: a.n]:
    with e.begin() as t:
        new = from_bytes(bytes(t.get(k)))
    old = None
    mrow = manifest[k.decode()]
    skey = mrow[10].encode() if len(mrow) > 10 else k       # complexes by assembly: the entry holds the assembly mrow[11]
    for s in src:
        with s.begin() as t:
            raw = t.get(skey)
        if raw is not None:
            old = from_bytes(bytes(raw))
            break
    assert old is not None
    if len(mrow) > 10:
        old = {mrow[11]: old[mrow[11]]}
    assert set(new) == set(old)
    for ak in new:
        nd, od = new[ak]["cifmol_attached_dict"], old[ak]["cifmol_attached_dict"]
        for part in ("atoms", "residues", "chains", "index_table"):
            assert json.dumps(nd[part], default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)) == json.dumps(od[part], default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)), (k, part)
        md = nd["metadata"]
        assert {k2: v for k2, v in md.items() if k2 in od["metadata"]} == od["metadata"], k
        cm = CIFMolAttached.from_dict(nd)
        _, s = secondary_structure(cm)
        m = manifest[k.decode()]
        ok = (abs(md["ss_helix"] - s["helix"]) < 1e-3 and abs(md["ss_strand"] - s["strand"]) < 1e-3 and abs(md["ss_coil"] - s["coil"]) < 1e-3
              and md["ss_max_coil_run"] == s["max_coil_run"] and md["length"] == s["n_assessed"] and md["fs_cluster"] == m[7] and md["fs_cluster_size"] == int(m[8]))
        if not ok:
            bad += 1
            print("mismatch", k, md, s)
print(f"checked {min(a.n, len(keys))} of {len(keys)} entries, mismatches {bad}")
print("metadata keys added:", sorted(set(md) - set(od["metadata"])))
