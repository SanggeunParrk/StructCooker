"""Find units whose polymer node holds residues of more than one chain (every unit, no sampling).

usage: merged_node_scan.py <db_root> <out.tsv> [workers=32]

A node's residues are stored chain by chain, so a protein / NA node whose polymer residues (parent letter
other than ``X``: a merged ligand has its own ``cif_idx``, and an unknown residue is ``X`` too, which makes this a
lower bound) have a ``cif_idx`` that is not strictly increasing, or whose largest ``cif_idx`` exceeds the key
sequence (+1), contains more than one chain: the node key is the key of one chain only. Cause: a ligand covalently bonded to two polymer chains used to pull
both polymer nodes into one (``adapt.build_assembly``, covalent merge). Output: one row per affected unit.
Run through sbatch (``merged_node_scan.sbatch``); reads every record of the DB.
"""

import csv
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from structcooker.instructions.unit.db import UnitDB

ROOT = sys.argv[1]
_db = None


def check(uid: str):
    global _db
    if _db is None:
        _db = UnitDB(ROOT)
    u = _db.get(uid)
    nodes, rs = u["nodes"], u["residues"]
    bad = []
    for i, kind in enumerate(nodes["kind"]):
        if kind not in ("protein", "na"):
            continue
        m = rs["node"] == i
        # a ligand merged into the node (covalent bond) carries its own cif_idx (1..); count polymer residues only
        cif = rs["cif_idx"][m & (rs["parent"] != "X")]
        seq_len = len(nodes["key"][i]) - 2
        if len(cif) and (np.any(np.diff(cif) <= 0) or cif.max() > seq_len + 1):
            bad.append((i, seq_len, len(cif), int(cif.max())))
    return uid, u["meta"]["pdb_id"], bad


def main() -> None:
    out = Path(sys.argv[2])
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 32
    with open(Path(ROOT) / "manifest.tsv") as f:
        ids = [r["unit_id"] for r in csv.DictReader(f, delimiter="\t")]
    n_bad = 0
    entries = set()
    with Pool(workers) as p, open(out, "w") as fo:
        fo.write("unit_id\tpdb_id\tnode\tkey_len\tn_res\tmax_cif_idx\n")
        for k, (uid, pdb, bad) in enumerate(
            p.imap_unordered(check, ids, chunksize=200), 1
        ):
            if bad:
                n_bad += 1
                entries.add(pdb)
                for node, kl, nr, mx in bad:
                    fo.write(f"{uid}\t{pdb}\t{node}\t{kl}\t{nr}\t{mx}\n")
            if k % 100000 == 0:
                print(f"{k}/{len(ids)} units, {n_bad} affected", flush=True)
    print(
        f"scanned {len(ids)} units: {n_bad} affected units in {len(entries)} entries -> {out}"
    )


if __name__ == "__main__":
    main()
