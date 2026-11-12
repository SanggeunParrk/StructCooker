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
    if a.pdb:  # a chosen entry (any shard), optionally only its first records: a scale and memory test
        wl = build_worklists(index, examples, {a.pdb})
        if a.max_records:
            keep = sorted(wl[a.pdb])[: a.max_records]
            wl[a.pdb] = {k: wl[a.pdb][k] for k in keep}
    pick = [a.pdb] if a.pdb else sorted(wl, key=lambda p: len(wl[p]))[: a.entries]
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


def _sweep_one(args: tuple[dict, str]) -> tuple[str, str | None, float]:
    rec, blob_dir = args
    t = time.time()
    try:
        value = _fresh(rec, Path(blob_dir))
        if not value["assembly_dict"]:
            return rec["key"], "empty assembly_dict", time.time() - t
    except Exception as e:  # noqa: BLE001 - a failing record is the finding
        return (
            rec["key"],
            f"{type(e).__name__}: {str(e).split(': ', 1)[-1]}"[:200],
            time.time() - t,
        )
    return rec["key"], None, time.time() - t


def cmd_sweep(a: argparse.Namespace) -> None:
    """Convert one random record of N random (entry, site, glycan) groups and group the failures by cause."""
    import multiprocessing as mp

    index, examples = _tables(a)
    on_disk = {p.name for p in Path(a.blob_dir).glob("*.bin")}
    rows = [
        (p, sh)
        for p, sh in zip(
            index.column("path").to_pylist(),
            index.column("shard").to_pylist(),
            strict=True,
        )
        if sh in on_disk
    ]
    groups: dict[str, list[str]] = {}
    for p, _ in rows:
        key = p.rsplit("/", 1)[-1].removesuffix(".cif")
        ref = key.rsplit("_conf", 1)[0]
        groups.setdefault(ref, []).append(key)
    rng = random.Random(a.seed)
    chosen = {
        rng.choice(v)
        for v in rng.sample(sorted(groups.values(), key=len), min(a.n, len(groups)))
    }
    pdbs = {k.split("_", 1)[0] for k in chosen}
    wl = build_worklists(index, examples, pdbs)
    recs = [r for e in wl.values() for r in e.values() if r["key"] in chosen]
    print(
        f"groups available {len(groups)}, converting {len(recs)} records of {len(pdbs)} entries",
        flush=True,
    )
    fail: dict[str, list[str]] = {}
    times = []
    with mp.get_context("fork").Pool(a.workers) as pool:
        for i, (key, err, dt) in enumerate(
            pool.imap_unordered(_sweep_one, [(r, str(a.blob_dir)) for r in recs])
        ):
            times.append(dt)
            if err:
                fail.setdefault(err, []).append(key)
            if (i + 1) % 200 == 0:
                print(
                    f"  {i + 1}/{len(recs)} failures {sum(map(len, fail.values()))}",
                    flush=True,
                )
    print(
        f"converted {len(recs) - sum(map(len, fail.values()))}/{len(recs)}; mean {statistics.mean(times):.2f} s/record"
    )
    for err, keys in sorted(fail.items(), key=lambda kv: -len(kv[1])):
        print(f"  {len(keys):5d}  {err}  e.g. {keys[:2]}")
    if fail:
        raise SystemExit(1)


def _decode_chunk(args: tuple[str, list[str]]) -> dict:
    """Decode every record of a key-contiguous chunk with one warm reader; count plain / delta values and bytes."""
    path, keys = args
    env = lmdb.open(path, readonly=True, lock=False, readahead=False)
    reader = DeltaReader(env)
    out = {
        "ok": 0,
        "bad": [],
        "plain": 0,
        "delta": 0,
        "plain_bytes": 0,
        "delta_bytes": 0,
        "secs": 0.0,
    }
    with env.begin() as txn:
        sizes = {k: len(bytes(txn.get(k.encode()))) for k in keys}
        kinds = {k: is_delta(bytes(txn.get(k.encode()))) for k in keys}
    for k in keys:
        t = time.time()
        try:
            v = reader.read(k)
            good = (
                bool(v.get("assembly_dict")) and "deposition_date" in v["metadata_dict"]
            )
        except Exception:  # noqa: BLE001 - a bad record is a finding
            good = False
        out["secs"] += time.time() - t
        if good:
            out["ok"] += 1
        elif len(out["bad"]) < 5:
            out["bad"].append(k)
        kind = "delta" if kinds[k] else "plain"
        out[kind] += 1
        out[f"{kind}_bytes"] += sizes[k]
    env.close()
    return out


def _equal_one(args: tuple[str, dict, str]) -> tuple[str, bool, float, bool]:
    """One record read with a cold reader against a fresh conversion of its CIF."""
    path, rec, blob_dir = args
    env = lmdb.open(path, readonly=True, lock=False, readahead=False)
    t = time.time()
    got = DeltaReader(env, cache_bytes=0).read(rec["key"])
    cold = time.time() - t
    with env.begin() as txn:
        delta = is_delta(bytes(txn.get(rec["key"].encode())))
    env.close()
    want = _fresh(rec, Path(blob_dir))
    return rec["key"], canonical_of(got) == canonical_of(want), cold, delta


def cmd_verify(a: argparse.Namespace) -> None:
    """Decode every record of a built DB, compare a random sample with fresh conversions, report sizes and timings."""
    import multiprocessing as mp

    index_path = Path(a.lmdb + ".index.tsv")
    keys = [
        ln.split("\t", 1)[0]
        for ln in index_path.read_text().splitlines()
        if not ln.startswith(("#", "key\t"))
    ]
    print(f"records in the index: {len(keys)}", flush=True)
    # 1. decode every record, in key-contiguous chunks (the records of an entry together: the entry base is decoded once)
    keys_sorted = sorted(keys)
    chunks, cur, last = [], [], None
    for k in keys_sorted:
        pdb = k.split("_", 1)[0]
        if len(cur) >= a.chunk and pdb != last:
            chunks.append(cur)
            cur = []
        cur.append(k)
        last = pdb
    chunks.append(cur)
    t0 = time.time()
    tot = {
        "ok": 0,
        "bad": [],
        "plain": 0,
        "delta": 0,
        "plain_bytes": 0,
        "delta_bytes": 0,
        "secs": 0.0,
    }
    with mp.get_context("fork").Pool(a.workers) as pool:
        for i, r in enumerate(
            pool.imap_unordered(_decode_chunk, [(a.lmdb, c) for c in chunks])
        ):
            for f in ("ok", "plain", "delta", "plain_bytes", "delta_bytes", "secs"):
                tot[f] += r[f]
            tot["bad"] += r["bad"]
            if (i + 1) % 20 == 0:
                print(
                    f"  chunks {i + 1}/{len(chunks)}  decoded {tot['ok']}  bad {len(tot['bad'])}  {time.time() - t0:.0f}s",
                    flush=True,
                )
    n = tot["ok"] + len(tot["bad"])
    print(
        f"decoded every record: {tot['ok']}/{len(keys)} ok (bad: {tot['bad'][:5]}) in {time.time() - t0:.0f}s wall"
    )
    print(
        f"  plain values {tot['plain']} ({tot['plain_bytes'] / 1e9:.2f} GB), delta values {tot['delta']} ({tot['delta_bytes'] / 1e9:.2f} GB)"
    )
    print(
        f"  warm decode (one reader per chunk): {1e3 * tot['secs'] / max(1, n):.1f} ms/record mean over {n} records"
    )
    # 2. a random sample against fresh conversions, each read with a cold reader
    index, examples = _tables(a)
    rng = random.Random(a.seed)
    chosen = set(rng.sample(keys, min(a.n, len(keys))))
    wl = build_worklists(index, examples, {k.split("_", 1)[0] for k in chosen})
    recs = [r for e in wl.values() for r in e.values() if r["key"] in chosen]
    bad, cold, n_delta = [], [], 0
    t1 = time.time()
    with mp.get_context("fork").Pool(a.workers) as pool:
        for i, (key, ok, c, delta) in enumerate(
            pool.imap_unordered(
                _equal_one, [(a.lmdb, r, str(a.blob_dir)) for r in recs]
            )
        ):
            cold.append(c)
            n_delta += delta
            if not ok:
                bad.append(key)
            if (i + 1) % 500 == 0:
                print(
                    f"  compared {i + 1}/{len(recs)}  mismatches {len(bad)}  {time.time() - t1:.0f}s",
                    flush=True,
                )
    cold.sort()
    print(
        f"equal to a fresh conversion: {len(recs) - len(bad)}/{len(recs)} ({n_delta} of them delta values); mismatches {bad[:5]}\n"
        f"cold decode ms/record (new reader, no cache): mean {1e3 * statistics.mean(cold):.1f}  median {1e3 * cold[len(cold) // 2]:.1f}  "
        f"p95 {1e3 * cold[int(0.95 * len(cold))]:.1f}  max {1e3 * cold[-1]:.1f}",
        flush=True,
    )
    if bad or tot["bad"] or tot["ok"] != len(keys):
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
    s.add_argument("--pdb", help="test this entry instead of the smallest of the shard")
    s.add_argument(
        "--max-records",
        type=int,
        default=0,
        help="with --pdb: only its first N records",
    )
    s.add_argument("--verify-n", type=int, default=40)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_sample)
    w = sp.add_parser("sweep")
    w.add_argument(
        "--n",
        type=int,
        default=3000,
        help="number of (entry, site, glycan) groups to sample",
    )
    w.add_argument("--workers", type=int, default=16)
    w.add_argument("--seed", type=int, default=0)
    w.set_defaults(fn=cmd_sweep)
    v = sp.add_parser("verify")
    v.add_argument("--lmdb", required=True)
    v.add_argument("--n", type=int, default=5000)
    v.add_argument("--workers", type=int, default=64)
    v.add_argument(
        "--chunk",
        type=int,
        default=4000,
        help="records per decode chunk (cut at an entry boundary)",
    )
    v.add_argument("--seed", type=int, default=0)
    v.set_defaults(fn=cmd_verify)
    a = ap.parse_args()
    a.fn(a)
    print(json.dumps({"done": True}))


if __name__ == "__main__":
    main()
