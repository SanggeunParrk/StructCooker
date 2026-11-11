"""List the cif_pdb entries in which a covalent ligand bridges two polymer chains (exact, every entry).

usage: bridge_entries_scan.py <out.tsv> [workers=32]

Runs ``adapt.build_assembly`` on the assemblies the unit DB reads (``select_assembly_keys``) and reports those with
``stats["n_bridge_skipped"] > 0``: the entries whose units were built with fused polymer nodes before the fix
(a ligand bonded to two polymer chains pulled both into one node, keyed by one chain's sequence). It is the
exact version of ``merged_node_scan.py``, which only sees the fused nodes that left a trace in the stored records.
Run through sbatch (``bridge_entries_scan.sbatch``); reads all of cif_pdb.
"""

import os
import sys
from multiprocessing import Pool
from pathlib import Path
from typing import cast

from datacooker.lmdb.sharded import open_env

from structcooker.instructions.mols import CIFMol
from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.instructions.unit.adapt import build_assembly
from structcooker.instructions.unit.config import load_config
from structcooker.instructions.unit.units import select_assembly_keys

DB = (
    os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")
    + "/lmdb/pdb/cif/cif_pdb.lmdb"
)
_env = None
_cfg = None


def check(key: bytes):
    global _env, _cfg
    if _env is None:
        _env = open_env(DB, readonly=True, lock=False, readahead=False)
        _cfg = load_config()
    with _env.begin() as txn:
        raw = txn.get(key)
    entry = from_bytes(bytes(raw))
    md = entry["metadata_dict"]
    hits = []
    for ck in select_assembly_keys(entry["assembly_dict"]):
        aid, mid, alt = ck.split("_")
        item = dict(entry["assembly_dict"][ck])
        item["metadata"] = {**md, "assembly_id": aid, "model_id": mid, "alt_id": alt}
        asm = build_assembly(
            CIFMol.from_dict(cast("dict", item)), _cfg.excluded, _cfg.min_heavy
        )
        if asm.stats.get("n_bridge_skipped", 0):
            hits.append((ck, asm.stats["n_bridge_skipped"]))
    return key.decode(), hits


def main() -> None:
    out = Path(sys.argv[1])
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 32
    env = open_env(DB, readonly=True, lock=False, readahead=False)
    with env.begin() as txn:
        keys = [bytes(k) for k in txn.cursor().iternext(values=False)]
    n_hit = 0
    with Pool(workers) as p, open(out, "w") as fo:
        fo.write("pdb_id\tassembly_key\tn_bridge_skipped\n")
        for k, (pdb, hits) in enumerate(p.imap_unordered(check, keys, chunksize=20), 1):
            for ck, n in hits:
                fo.write(f"{pdb}\t{ck}\t{n}\n")
            n_hit += bool(hits)
            if k % 20000 == 0:
                print(
                    f"{k}/{len(keys)} entries, {n_hit} with a bridging ligand",
                    flush=True,
                )
    print(f"scanned {len(keys)} entries: {n_hit} with a bridging ligand -> {out}")


if __name__ == "__main__":
    main()
