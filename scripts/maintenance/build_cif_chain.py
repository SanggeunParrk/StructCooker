"""Pre-build a per-chain cif LMDB for fast template lookups.

Reads a cif LMDB ({assembly_dict, metadata_dict} per PDB id) and writes one
entry per chain, keyed ``<pdbid>_<chain>`` (lower-cased pdb id), value = the
chain's extracted BioMol dict (best-occupancy assembly, via
``cif_record_to_chains``). Template ingest can then resolve each hit with a
light keyed read instead of re-decoding + rebuilding every assembly per hit.

Sharded for SLURM arrays: ``--shard-idx i --n-shards N`` processes the i-th
1/N slice of source keys into ``<out>_shard<i>.lmdb`` (merge afterwards).

    pixi run python -m scripts.maintenance.build_cif_chain --cif-db CIF --out-db OUT
        --shard-idx 0 --n-shards 20 --n-jobs 16
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import lmdb
from joblib import Parallel, delayed

from structcooker.instructions.readers.io import load_bytes
from structcooker.instructions.transforms.codecs import to_bytes
from structcooker.instructions.transforms.template import cif_record_to_chains

_BATCH = 200  # PDB ids processed (in parallel) per write transaction
# OOM guard: skip records whose largest assembly exceeds this atom count
# (a handful of giant assemblies would otherwise OOM a worker). Override via env.
_MAX_ATOMS = int(os.environ.get("DC_CIFCHAIN_MAX_ATOMS", "1500000") or "1500000")


def _chains_for(pdb_id: str, raw: bytes) -> list[tuple[bytes, bytes]]:
    """Return [(key, value)] per chain for one cif record, or [] on failure/skip."""
    try:
        record = load_bytes(raw)
        if _MAX_ATOMS > 0:
            max_atoms = max(
                (len(a["atoms"]["nodes"]["id"]["value"]) for a in record["assembly_dict"].values()),
                default=0,
            )
            if max_atoms > _MAX_ATOMS:
                return []
        chains = cif_record_to_chains(record)
    except Exception:  # noqa: BLE001 - skip unparsable records, report via count
        return []
    return [
        (f"{pdb_id.lower()}_{base}".encode(), to_bytes(biomoldict))
        for base, biomoldict in chains.items()
    ]


def build(cif_db: Path, out_db: Path, shard_idx: int, n_shards: int, n_jobs: int) -> None:
    """Build (a shard of) the per-chain cif LMDB."""
    env_in = lmdb.open(str(cif_db), readonly=True, lock=False, subdir=True, max_dbs=0)
    out_db.parent.mkdir(parents=True, exist_ok=True)
    env_out = lmdb.open(str(out_db), map_size=int(3e11), subdir=True)

    pdbs = written = failed = 0
    batch: list[tuple[str, bytes]] = []

    def flush(items: list[tuple[str, bytes]]) -> None:
        nonlocal written, failed
        results = Parallel(n_jobs=n_jobs, prefer="processes")(
            delayed(_chains_for)(pid, raw) for pid, raw in items
        )
        with env_out.begin(write=True) as txn:
            for entries in results:
                if not entries:
                    failed += 1
                    continue
                for key, value in entries:
                    txn.put(key, value, overwrite=True)
                    written += 1

    with env_in.begin() as txn:
        # Keys only first (cheap), shard, then fetch values for our slice only --
        # avoids every shard reading every (large) record in the DB.
        all_keys = list(txn.cursor().iternext(keys=True, values=False))
        shard_keys = [bytes(all_keys[i]) for i in range(shard_idx, len(all_keys), n_shards)]
        for key in shard_keys:
            batch.append((key.decode(), bytes(txn.get(key))))
            pdbs += 1
            if len(batch) >= _BATCH:
                flush(batch)
                batch = []
        if batch:
            flush(batch)

    env_in.close()
    env_out.close()
    sys.stdout.write(
        f"shard {shard_idx}/{n_shards}: pdbs={pdbs} chains_written={written} failed={failed}\n",
    )


def main() -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Build a per-chain cif LMDB.")
    parser.add_argument("--cif-db", type=Path, required=True)
    parser.add_argument("--out-db", type=Path, required=True)
    parser.add_argument("--shard-idx", type=int, default=0)
    parser.add_argument("--n-shards", type=int, default=1)
    parser.add_argument("--n-jobs", type=int, default=8)
    args = parser.parse_args()
    out = args.out_db
    if args.n_shards > 1:
        out = out.with_name(f"{out.stem}_shard{args.shard_idx}{out.suffix}")
    build(args.cif_db, out, args.shard_idx, args.n_shards, args.n_jobs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
