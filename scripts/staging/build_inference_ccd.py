"""Build the inference CCD (CCDMol graphs) from the release CCD (schema G).

    python scripts/staging/build_inference_ccd.py SRC_CCD.lmdb DST.lmdb [--jobs N]

SRC is the release ``ccd/biomol_CCD.lmdb``; DST (new, never replaced) gets one CCDMol per component
(:mod:`structcooker.instructions.transforms.ccdmol`). Components the CCD defines without atoms (UNL)
have no graph and are reported, not written. Run as a batch job.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import lmdb

from structcooker.instructions.transforms.ccdmol import convert_chem_comp_to_ccdmol
from structcooker.instructions.transforms.codecs import from_bytes, to_bytes


def _convert(args: tuple[str, list[bytes]]) -> list[tuple[bytes, bytes | None, str | None]]:
    src, keys = args
    out = []
    with lmdb.open(src, readonly=True, lock=False) as env, env.begin() as txn:
        for key in keys:
            try:
                out.append((key, to_bytes(convert_chem_comp_to_ccdmol(from_bytes(txn.get(key)), key.decode())), None))
            except Exception as exc:  # noqa: BLE001 - reported per component
                out.append((key, None, f"{type(exc).__name__}: {exc}"[:160]))
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("src", type=Path)
    parser.add_argument("dst", type=Path)
    parser.add_argument("--jobs", type=int, default=32)
    args = parser.parse_args()
    if args.dst.exists():
        raise SystemExit(f"{args.dst} exists; build into a new path")
    with lmdb.open(str(args.src), readonly=True, lock=False) as env, env.begin() as txn:
        keys = [bytes(k) for k in txn.cursor().iternext(values=False)]
    shards = [(str(args.src), keys[i::args.jobs]) for i in range(args.jobs)]
    with ProcessPoolExecutor(args.jobs) as pool:
        results = [r for part in pool.map(_convert, shards) for r in part]
    failed = [(k.decode(), err) for k, v, err in results if v is None]
    with lmdb.open(str(args.dst), map_size=100 * 1024**3) as env, env.begin(write=True) as txn:
        for key, value, _ in sorted(results):
            if value is not None:
                txn.put(key, value)
        written = txn.stat()["entries"]
    print(f"components {len(keys)}, written {written}, not converted {len(failed)}: {failed[:10]}")


if __name__ == "__main__":
    main()
