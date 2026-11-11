#!/usr/bin/env python3
"""Check the Re-Glyco delta DB (``glycan/reglyco``) against a fresh conversion of the downloaded CIFs.

  sample   encode a few whole entries of the downloaded component into a scratch LMDB with the stage's own code
           (``readers.reglyco.explode_reglyco_entry``), then check the records like ``verify`` does. A dry run of the
           stage on a handful of entries.
  verify   read N random records of a built DB through ``DeltaReader`` (cold: a new reader per record; warm: the records
           of one entry in a row), convert the same record again from its CIF and require the two values to be equal
           (canonical serialised bytes). Also: every key decodes, a plain ``from_bytes`` refuses a delta value, and the
           sizes and decode times are reported.

Heavy (it converts CIFs, 1-2 s each): run through sbatch (``reglyco_check.sbatch``).
"""

from __future__ import annotations

import argparse
import json
import os
import random
import statistics
import time
from pathlib import Path

import lmdb
import pyarrow.parquet as pq

from structcooker.instructions.readers.reglyco import (
    BLOB_ENV,
    convert_record,
    explode_reglyco_entry,
    read_blob,
)
from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.instructions.transforms.delta_codec import (
    DeltaReader,
    canonical_of,
    is_delta,
)
from structcooker.instructions.transforms.reglyco_worklist import build_worklists

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
COMP = (
    ROOT
    / "materials/raw/rfd4_proteina/train_datasets/streaming-latest/synthetic/glycans/2026_09_27_reglyco"
)


def _tables(a: argparse.Namespace):
    index = pq.read_table(
        a.index, columns=["path", "shard", "offset", "length", "cif_sha256"]
    )
    examples = pq.read_table(a.examples)
    return index, examples


def _fresh(wl_rec: dict, blob_dir: Path) -> dict:
    meta = {
        "deposition_date": wl_rec["deposition_date"],
        "resolution": wl_rec.get("resolution"),
    }
    return convert_record(read_blob(blob_dir, wl_rec), wl_rec["key"], meta)


def cmd_sample(a: argparse.Namespace) -> None:
    index, examples = _tables(a)
    shard_rows = [
        i for i, s in enumerate(index.column("shard").to_pylist()) if s == a.shard
    ]
    pdbs = sorted(
        {index.column("path")[i].as_py().rsplit("/", 1)[-1][:4] for i in shard_rows}
    )
    wl = build_worklists(index, examples, set(pdbs))
    pick = sorted(wl, key=lambda p: len(wl[p]))[: a.entries]
    print(
        f"entries of {a.shard}: {len(pdbs)}; testing {pick} with {[len(wl[p]) for p in pick]} records",
        flush=True,
    )
    os.environ[BLOB_ENV] = str(a.blob_dir)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(out / "scratch.lmdb"), map_size=1 << 33)
    t0 = time.time()
    n = total = 0
    for pdb in pick:
        t1 = time.time()
        sub = explode_reglyco_entry({"records": wl[pdb]})
        with env.begin(write=True) as txn:
            for k, v in sub.items():
                value = v["reglyco_record"]
                txn.put(f"{pdb}_{k}".encode(), value)
                total += len(value)
                n += 1
        print(f"  {pdb}: {len(sub)} records in {time.time() - t1:.0f}s", flush=True)
    print(
        f"encoded {n} records: {total / 1e6:.2f} MB ({total / n / 1e3:.1f} KB/record), {(time.time() - t0) / n:.2f} s/record"
    )
    env.close()
    a.lmdb, a.n = str(out / "scratch.lmdb"), a.verify_n
    a.wl = {k: r for p in pick for k, r in ((x["key"], x) for x in wl[p].values())}
    _verify(a, a.wl)


def _verify(a: argparse.Namespace, wl_by_key: dict[str, dict]) -> None:
    env = lmdb.open(a.lmdb, readonly=True, lock=False)
    with env.begin() as txn:
        keys = [bytes(k).decode() for k, _ in txn.cursor()]
        sizes = [len(bytes(v)) for _, v in txn.cursor()]
    n_delta = sum(1 for _, v in env.begin().cursor() if is_delta(bytes(v)))
    print(
        f"records {len(keys)}  stored {sum(sizes) / 1e6:.2f} MB  deltas {n_delta}  plain {len(keys) - n_delta}"
    )
    rng = random.Random(0)
    chosen = rng.sample(keys, min(a.n, len(keys)))
    blob_dir = Path(os.environ[BLOB_ENV])
    # 1. a plain from_bytes refuses a delta value by name
    with env.begin() as txn:
        deltas = [k for k in chosen if is_delta(bytes(txn.get(k.encode())))]
        plains = [k for k in chosen if k not in set(deltas)]
        refused = 0
        for k in deltas[:20]:
            try:
                from_bytes(bytes(txn.get(k.encode())))
            except ValueError as e:
                refused += "read_record" in str(e)
        print(
            f"plain from_bytes refuses a delta by name: {refused}/{min(20, len(deltas))}"
        )
        for k in plains[:5]:
            from_bytes(
                bytes(txn.get(k.encode()))
            )  # an entry base reads as an ordinary value
    # 2. decode == fresh conversion; cold and warm decode times
    cold, warm, bad = [], [], []
    warm_reader = DeltaReader(env)
    for k in sorted(chosen):  # sorted: the records of an entry come in a row
        t = time.time()
        got = DeltaReader(env, cache_bytes=0).read(k)
        cold.append(time.time() - t)
        t = time.time()
        got_w = warm_reader.read(k)
        warm.append(time.time() - t)
        want = _fresh(wl_by_key[k], blob_dir)
        if canonical_of(got) != canonical_of(want) or canonical_of(
            got_w
        ) != canonical_of(want):
            bad.append(k)
    print(
        f"equal to a fresh conversion: {len(chosen) - len(bad)}/{len(chosen)}; first mismatches {bad[:3]}\n"
        f"decode ms/record: cold mean {1e3 * statistics.mean(cold):.1f} (p95 {1e3 * sorted(cold)[int(0.95 * len(cold))]:.1f}), "
        f"warm mean {1e3 * statistics.mean(warm):.1f} (p95 {1e3 * sorted(warm)[int(0.95 * len(warm))]:.1f}); "
        f"base decodes {warm_reader.base_decodes}"
    )
    if bad:
        raise SystemExit(1)


def cmd_verify(a: argparse.Namespace) -> None:
    index, examples = _tables(a)
    with lmdb.open(a.lmdb, readonly=True, lock=False).begin() as txn:
        keys = [bytes(k).decode() for k, _ in txn.cursor()]
    rng = random.Random(0)
    chosen = set(rng.sample(keys, min(a.n, len(keys))))
    pdbs = {k.split("_", 1)[0] for k in chosen}
    wl = build_worklists(index, examples, pdbs)
    by_key = {r["key"]: r for e in wl.values() for r in e.values()}
    os.environ[BLOB_ENV] = str(a.blob_dir)
    # decode every key once (all of them must decode), then compare the sample with fresh conversions
    env = lmdb.open(a.lmdb, readonly=True, lock=False)
    reader = DeltaReader(env)
    t = time.time()
    fail = []
    for i, k in enumerate(keys):
        try:
            v = reader.read(k)
            if "assembly_dict" not in v or not v["assembly_dict"]:
                fail.append(k)
        except Exception:  # noqa: BLE001 - a bad record is a finding
            fail.append(k)
        if (i + 1) % 100000 == 0:
            print(
                f"  decoded {i + 1}/{len(keys)}  {time.time() - t:.0f}s  failures {len(fail)}",
                flush=True,
            )
    print(
        f"decoded every record: {len(keys) - len(fail)}/{len(keys)} ok in {time.time() - t:.0f}s"
    )
    env.close()
    a.n = len(chosen)
    _verify(
        argparse.Namespace(lmdb=a.lmdb, n=len(chosen)), {k: by_key[k] for k in chosen}
    )
    if fail:
        print("failures:", fail[:5])
        raise SystemExit(1)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--index", default=str(COMP / "blob/index.parquet"))
    ap.add_argument("--examples", default=str(COMP / "examples.parquet"))
    ap.add_argument("--blob-dir", type=Path, default=COMP / "blob/data")
    sp = ap.add_subparsers(required=True)
    s = sp.add_parser("sample")
    s.add_argument("--shard", default="rs_00005.bin")
    s.add_argument("--entries", type=int, default=1)
    s.add_argument("--verify-n", type=int, default=40)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_sample)
    v = sp.add_parser("verify")
    v.add_argument("--lmdb", required=True)
    v.add_argument("--n", type=int, default=5000)
    v.set_defaults(fn=cmd_verify)
    a = ap.parse_args()
    a.fn(a)
    print(json.dumps({"done": True}))


if __name__ == "__main__":
    main()
