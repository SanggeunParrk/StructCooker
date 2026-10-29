"""Fragment template DB (k-mer pairs with coordinates) from the training set (sbatch only).

usage: build_fragdb.py --out DIR [--k 7] [--workers 56] [--limit N] [--seed 0] [--flush 2000000]
Source: $BIOMOL_ROOT/lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb (first assembly of each entry).
Each worker handles a fixed share of the entries and writes parts DIR/part_<worker>_<n>/<field>.npy; see
structcooker.fragdb.extract for the content of a pair. ``pdb`` is the 4-character id packed into a uint32
(little endian: ``int.to_bytes(4, 'little').decode()``), ``assembly`` its number, ``date`` days since 1970-01-01. Release date and resolution come from metadata/cif_metadata.tsv.
"""
import argparse
import datetime
import os

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"

import csv  # noqa: E402
import multiprocessing as mp  # noqa: E402
import random  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import lmdb  # noqa: E402
import numpy as np  # noqa: E402

from structcooker.fragdb.extract import extract_pairs  # noqa: E402
from structcooker.instructions.transforms.codecs import from_bytes  # noqa: E402
from structcooker.mols import CIFMolAttached  # noqa: E402

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
DB = ROOT / "lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb"
EPOCH = datetime.date(1970, 1, 1)
META: dict[str, tuple[int, float]] = {}


def pick(keys):
    def k(c):
        a, m, alt = c.split("_")
        return (int(a), int(m), alt != ".", alt)
    return min(keys, key=k)


def flush(out, sid, n, bufs):
    d = Path(out) / f"part_{sid:03d}_{n:03d}"
    d.mkdir(parents=True, exist_ok=True)
    for f, arrs in bufs.items():
        np.save(d / f"{f}.npy", np.concatenate(arrs))
    return sum(len(a) for a in bufs["inter"])


def run_shard(task):
    sid, keys, out, k, flush_at = task
    env = lmdb.open(str(DB), readonly=True, lock=False, readahead=False, max_readers=512)
    bufs = {f: [] for f in ("seq", "coords", "center", "chain", "cluster", "inter", "n_atom_pairs", "pdb", "assembly", "date", "resolution")}
    n_buf = n_pairs = n_entries = n_fail = part = 0
    fails = []
    t0 = time.time()
    for idx, key in enumerate(keys):
        try:
            with env.begin() as txn:
                raw = txn.get(key.encode())
            rec = from_bytes(bytes(raw))
            ck = pick(rec)
            cm = CIFMolAttached.from_dict(rec[ck]["cifmol_attached_dict"])
            res = extract_pairs(cm, k=k)
            if res is not None:
                n = len(res["inter"])
                date, resolution = META.get(key, (0, np.nan))
                for f in res:
                    bufs[f].append(res[f])
                bufs["pdb"].append(np.full(n, int.from_bytes(key.encode()[:4].ljust(4), "little"), dtype=np.uint32))   # 4 characters packed
                bufs["assembly"].append(np.full(n, int(ck.split("_")[0]) if ck.split("_")[0].isdigit() else 0, dtype=np.uint16))
                bufs["date"].append(np.full(n, date, dtype=np.int32))
                bufs["resolution"].append(np.full(n, resolution, dtype=np.float16))
                n_buf += n; n_pairs += n; n_entries += 1
                if n_buf >= flush_at:
                    flush(out, sid, part, bufs); part += 1; n_buf = 0
                    bufs = {f: [] for f in bufs}
        except Exception:  # noqa: BLE001
            n_fail += 1
            fails.append((key, traceback.format_exc().strip().splitlines()[-1][:160]))
        if (idx + 1) % 1000 == 0:
            print(f"shard {sid}: {idx + 1}/{len(keys)} entries, {n_pairs} pairs, {time.time() - t0:.0f}s", flush=True)
    if n_buf:
        flush(out, sid, part, bufs)
    return sid, len(keys), n_entries, n_pairs, n_fail, fails[:3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=7)
    ap.add_argument("--workers", type=int, default=56)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--flush", type=int, default=2_000_000)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for r in csv.DictReader(open(ROOT / "metadata/cif_metadata.tsv"), delimiter="\t"):
        pid = r["cif_id"].split("_")[0]
        if pid not in META and r["release_date"]:
            META[pid] = ((datetime.date.fromisoformat(r["release_date"]) - EPOCH).days, float(r["resolution"]) if r["resolution"] not in ("", "None") else np.nan)
    sizes = {}
    for line in open(str(DB) + ".index.tsv"):
        if line.startswith(("#", "key")):
            continue
        key, b = line.split("\t")[:2]
        sizes[key] = int(b)
    keys = list(sizes)
    random.Random(a.seed).shuffle(keys)
    if a.limit:
        keys = keys[: a.limit]
    keys.sort(key=lambda x: -sizes[x])
    tasks = [(s, keys[s :: a.workers], str(out), a.k, a.flush) for s in range(a.workers)]
    print(f"{len(keys)} entries, {a.workers} shards, k={a.k}", flush=True)
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as pool:
        res = list(pool.imap_unordered(run_shard, tasks))
    n_pairs = sum(r[3] for r in res)
    print(f"done in {time.time() - t0:.0f}s: {sum(r[2] for r in res)} of {sum(r[1] for r in res)} entries gave pairs, {n_pairs} pairs, "
          f"{sum(r[4] for r in res)} entries failed", flush=True)
    for r in res:
        for f in r[5]:
            print("  fail:", f, flush=True)
    (out / "SUMMARY.txt").write_text(f"k={a.k}\nentries={sum(r[1] for r in res)}\nentries_with_pairs={sum(r[2] for r in res)}\npairs={n_pairs}\n")


if __name__ == "__main__":
    main()
