"""Check a built design DB against its source and against a fresh statistics computation (sbatch only).

usage: validate_db.py DESIGN_DB SOURCE_DB [--n 500]
For a sample of the design DB: the structure equals the source record's, the metadata keeps the source's and adds the
design keys, and the statistics in the metadata equal ``instructions.design.stats.item_stats`` of the structure.
"""

import argparse
import random

from datacooker.lmdb.sharded import open_env

from structcooker.instructions.design.stats import item_stats, load_item
from structcooker.instructions.transforms.codecs import from_bytes

ap = argparse.ArgumentParser()
ap.add_argument("db")
ap.add_argument("src")
ap.add_argument("--n", type=int, default=500)
a = ap.parse_args()
kw = {"readonly": True, "lock": False}
bad = 0
with (
    open_env(a.db, **kw) as env,
    open_env(a.src, **kw) as src,
    env.begin() as t,
    src.begin() as s,
):
    keys = list(t.cursor().iternext(values=False))
    random.Random(0).shuffle(keys)
    for k in keys[: a.n]:
        raw = bytes(t.get(k))
        new, old = from_bytes(raw), from_bytes(bytes(s.get(k)))
        assert set(new) == set(old), k
        for ak in new:
            nd, od = new[ak]["cifmol_attached_dict"], old[ak]["cifmol_attached_dict"]
            for part in ("atoms", "residues", "chains", "index_table"):
                assert repr(nd[part]) == repr(od[part]), (k, part)
            md = nd["metadata"]
            assert {x: y for x, y in md.items() if x in od["metadata"]} == od[
                "metadata"
            ], k
        cols, _ = item_stats(load_item(raw))
        md = new[sorted(new)[0]]["cifmol_attached_dict"]["metadata"]
        got = [
            md["length"],
            md["mean_plddt"],
            md["ss_helix"],
            md["ss_strand"],
            md["ss_coil"],
            md["ss_max_coil_run"],
        ]
        if any(abs(float(x) - float(y)) > 1e-3 for x, y in zip(got, cols, strict=True)):
            bad += 1
            print("mismatch", k, got, cols)
print(f"checked {min(a.n, len(keys))} of {len(keys)} entries, mismatches {bad}")
