"""Compare structcooker.design.ss with mkdssp on a sample of PDB items (sbatch only).

usage: check_ss.py --n 100 [--seed 0] [--work DIR]
Writes N/CA/C/O of every protein residue as a PDB file, runs ``mkdssp`` (legacy output) and compares the three-state codes.
"""
import argparse
import os
import random
import subprocess
import tempfile
from pathlib import Path

import lmdb
import numpy as np

from structcooker.design.export import write_backbone_pdb
from structcooker.design.ss import secondary_structure
from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.mols import CIFMolAttached

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
DB = ROOT / "lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb"
AA = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
MAP = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E"}


def write_pdb(cm, path):
    return write_backbone_pdb(cm, path)


def read_dssp(path):
    out, on = [], False
    for ln in Path(path).read_text().splitlines():
        if ln.startswith("  #  RESIDUE"):
            on = True
        elif on and len(ln) > 16 and ln[13] != "!":
            out.append(MAP.get(ln[16], "C"))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--work", default=None)
    a = ap.parse_args()
    work = Path(a.work or tempfile.mkdtemp())
    work.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(DB), readonly=True, lock=False, max_readers=8)
    with env.begin() as txn:
        keys = [k for k, _ in zip(txn.cursor().iternext(values=False), range(20000))]
    random.Random(a.seed).shuffle(keys)
    agree = tot = done = 0
    fr_err, rows = [], []
    for key in keys:
        if done >= a.n:
            break
        with env.begin() as txn:
            rec = from_bytes(txn.get(key))
        akey = sorted(rec)[0]
        cm = CIFMolAttached.from_dict(rec[akey]["cifmol_attached_dict"]) if hasattr(CIFMolAttached, "from_dict") else None
        if cm is None:
            raise SystemExit("cannot build CIFMolAttached from record; adapt loader")
        s3, summ = secondary_structure(cm)
        if summ["n_assessed"] < 30:
            continue
        pdb = work / f"{key.decode()}.pdb"
        use = write_pdb(cm, pdb)
        res = subprocess.run(["pixi", "run", "mkdssp", "--output-format", "dssp", str(pdb), str(work / "o.dssp")], capture_output=True, text=True)
        if res.returncode:
            print("mkdssp failed", key, res.stderr[-200:])
            continue
        ref = read_dssp(work / "o.dssp")
        mine = s3[use]
        if len(ref) != len(mine):
            print("length mismatch", key, len(ref), len(mine))
            continue
        agree += int((ref == mine).sum())
        tot += len(ref)
        f = [float((ref == c).mean()) for c in "HEC"]
        fr_err.append([summ["helix"] - f[0], summ["strand"] - f[1], summ["coil"] - f[2]])
        done += 1
    e = np.array(fr_err)
    print(f"items {done}  residues {tot}  per-residue 3-state agreement {agree / tot:.4f}")
    print("fraction error mean (mine - dssp) H/E/C", e.mean(0).round(4), " mean |err|", np.abs(e).mean(0).round(4), " max |err|", np.abs(e).max(0).round(3))


if __name__ == "__main__":
    main()
