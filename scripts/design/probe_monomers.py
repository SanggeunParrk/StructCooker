"""Small end-to-end test on the OpenFold short distillation set (sbatch only): per-item stats + foldseek clustering.

usage: probe_monomers.py --n 5000 --work DIR [--seed 0]
Per item: length, mean pLDDT (b-factor), secondary structure fractions. Then ``foldseek easy-cluster`` on the backbone PDBs.
"""
import argparse
import csv
import os
import random
import subprocess
import time
from pathlib import Path

import lmdb
import numpy as np

from structcooker.design.export import write_backbone_pdb
from structcooker.design.ss import secondary_structure
from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.mols import CIFMolAttached

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
DB = ROOT / "lmdb/distillation_short/cif_attached/cif_short_attached.lmdb"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=5000)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    work = Path(a.work)
    (work / "pdb").mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(DB), readonly=True, lock=False, subdir=True)
    with env.begin() as txn:
        keys = list(txn.cursor().iternext(values=False))
    random.Random(a.seed).shuffle(keys)
    t0 = time.time()
    rows = []
    with env.begin() as txn:
        for key in keys[: a.n]:
            rec = from_bytes(txn.get(key))
            cm = CIFMolAttached.from_dict(rec[sorted(rec)[0]]["cifmol_attached_dict"])
            _, s = secondary_structure(cm)
            name = key.decode()
            use = write_backbone_pdb(cm, work / "pdb" / f"{name}.pdb")
            plddt = float(np.mean(cm.atoms.b_factor.value))
            rows.append({"id": name, "length": s["n_assessed"], "plddt": round(plddt, 2), **{k: round(s[k], 4) for k in ("helix", "strand", "coil")}, "max_coil_run": s["max_coil_run"]})
    print(f"stats {len(rows)} items in {time.time() - t0:.0f}s", flush=True)
    with open(work / "stats.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    t0 = time.time()
    r = subprocess.run(["pixi", "run", "foldseek", "easy-cluster", str(work / "pdb"), str(work / "clu"), str(work / "tmp"), "--threads", "4",
                        "-c", "0.9", "--min-seq-id", "0.0", "--alignment-type", "1", "--tmscore-threshold", "0.5"],
                       capture_output=True, text=True)
    print("foldseek", r.returncode, f"{time.time() - t0:.0f}s", r.stderr[-300:], flush=True)
    cl = np.loadtxt(work / "clu_cluster.tsv", dtype=str, delimiter="\t")
    reps = len(set(cl[:, 0]))
    print(f"{len(cl)} members, {reps} clusters")


if __name__ == "__main__":
    main()
