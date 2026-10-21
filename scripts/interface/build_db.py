"""Build the PDB interface-unit DB: one LMDB shard per worker + manifest.tsv (sbatch only).

usage: build_db.py --out DIR [--workers 96] [--min-mb 0] [--max-mb 30] [--limit N] [--seed 0]
Reads $BIOMOL_ROOT/lmdb/pdb/cif/cif_pdb.lmdb (waters included). Resumable: finished entries are
listed in DIR/done.tsv and skipped on restart. Per entry: first model of every assembly (altloc '.'
preferred), units deduplicated across the entry's assemblies by (node sequences, contact size).
"""
import argparse
import os

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"

import atexit  # noqa: E402
import multiprocessing as mp  # noqa: E402
import random  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import lmdb  # noqa: E402
import numpy as np  # noqa: E402

from structcooker.instructions.transforms.codecs import (  # noqa: E402
    from_bytes,
    to_bytes,
)
from structcooker.interface.antibody import AntibodyClassifier  # noqa: E402
from structcooker.interface.config import load_config  # noqa: E402
from structcooker.interface.pipeline import extract_units  # noqa: E402
from structcooker.interface.record import (  # noqa: E402
    asm_fields,
    manifest_row,
    unit_to_record,
)
from structcooker.mols import CIFMol  # noqa: E402

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
CFG = AB = ENV_IN = ENV_OUT = None
WID = -1
SHARD_BASE = 0
MAP_SIZE = 64 << 30


def select_keys(keys):
    """First model per assembly; the no-altloc entry if present, else the first altloc."""
    best = {}
    for k in sorted(keys):
        a, m, alt = k.split("_")
        cur = best.get(a)
        if cur is None or int(m) < int(cur[1]) or (m == cur[1] and alt == "."):
            best[a] = (k, m)
    return [v[0] for _, v in sorted(best.items(), key=lambda kv: kv[0])]


def init_worker(counter, out_dir, ab_cache, shard_base):
    global CFG, AB, ENV_IN, ENV_OUT, WID, SHARD_BASE
    with counter.get_lock():
        WID = counter.value
        counter.value += 1
    SHARD_BASE = shard_base
    CFG = load_config()
    AB = AntibodyClassifier(ab_cache)  # cache misses (rare) fall back to a live ANARCI call
    ENV_IN = lmdb.open(str(ROOT / "lmdb/pdb/cif/cif_pdb.lmdb"), readonly=True, lock=False, readahead=False, max_readers=512)
    shard = Path(out_dir) / "shards" / f"shard_{SHARD_BASE + WID:03d}.lmdb"
    shard.parent.mkdir(parents=True, exist_ok=True)
    ENV_OUT = lmdb.open(str(shard), map_size=MAP_SIZE, subdir=True, metasync=False, sync=False, max_dbs=0)
    atexit.register(ENV_OUT.sync)


def process_entry(key):
    with ENV_IN.begin() as txn:
        raw = txn.get(key.encode())
    rec = from_bytes(bytes(raw))
    del raw
    md = rec["metadata_dict"]
    res = md.get("resolution")
    base_meta = {"pdb_id": key, "resolution": None if res is None or not np.isfinite(float(res)) else float(res),
                 "release_date": md.get("release_date"), "deposition_date": md.get("deposition_date")}
    out_rows, seen, n_skipped = [], set(), 0
    puts = []
    for ck in select_keys(rec["assembly_dict"]):
        aid, mid, alt = ck.split("_")
        item = dict(rec["assembly_dict"][ck])
        m = dict(md)
        m.update(assembly_id=aid, model_id=mid, alt_id=alt)
        item["metadata"] = m
        cm = CIFMol.from_dict(item)
        asm, units, duplex = extract_units(cm, CFG, AB)
        if sum(k in ("protein", "na") for k in asm.kinds) > CFG.unit.max_polymer_chains_in_assembly:
            n_skipped += 1
            del cm, asm
            continue
        if not units:
            del cm, asm
            continue
        f = asm_fields(cm, asm)
        idx = 0
        for u in units:
            sig = (tuple(sorted(asm.keys[x] for x in u.nodes)), round(u.weight / 10))
            if sig in seen:
                continue
            seen.add(sig)
            meta = dict(base_meta, assembly_id=aid, model_id=mid, alt_id=alt)
            r = unit_to_record(asm, f, u, duplex, CFG, meta, asm.stats["ab_class"], idx)
            idx += 1
            row = manifest_row(r)
            row["shard"] = SHARD_BASE + WID
            out_rows.append(row)
            puts.append((r["meta"]["unit_id"].encode(), to_bytes(r, level=5)))
        del cm, asm, f
    if puts:
        with ENV_OUT.begin(write=True) as txn:
            for k, v in puts:
                txn.put(k, v)
    return out_rows, n_skipped


def run_chunk(keys):
    rows, done, fails = [], [], []
    for k in keys:
        try:
            r, skipped = process_entry(k)
            rows += r
            done.append((k, len(r), skipped))
        except Exception:  # noqa: BLE001
            fails.append((k, traceback.format_exc().strip().splitlines()[-1][:300]))
    ENV_OUT.sync()
    return rows, done, fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=96)
    ap.add_argument("--min-mb", type=float, default=0)
    ap.add_argument("--max-mb", type=float, default=30)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=8)
    ap.add_argument("--antibody-cache", default=None)
    ap.add_argument("--keys", default=None, help="file with one PDB key per line (overrides the size/limit selection)")
    a = ap.parse_args()
    out = Path(a.out)
    if (out / "interface_pdb.lmdb" / "shards.json").exists():
        raise SystemExit(f"{out} is already finalized (interface_pdb.lmdb published); build into a new directory")
    out.mkdir(parents=True, exist_ok=True)
    index = ROOT / "lmdb/pdb/cif/cif_pdb.lmdb.index.tsv"
    sizes = {}
    for line in open(index):
        if line.startswith(("#", "key")):
            continue
        k, b = line.split("\t")[:2]
        sizes[k] = int(b)
    done_p = out / "done.tsv"
    done = set(l.split("\t")[0] for l in open(done_p)) if done_p.exists() else set()
    if a.keys:
        keys = [l.strip() for l in open(a.keys) if l.strip()]
    else:
        keys = [k for k, b in sizes.items() if a.min_mb * 1e6 <= b < a.max_mb * 1e6]
        random.Random(a.seed).shuffle(keys)
        if a.limit:
            keys = keys[: a.limit]
    keys = [k for k in keys if k not in done]
    keys.sort(key=lambda k: -sizes[k])  # heaviest first: balances the tail
    chunks = [keys[i : i + a.chunk] for i in range(0, len(keys), a.chunk)]
    shard_base = len(list((out / "shards").glob("shard_*.lmdb"))) if (out / "shards").exists() else 0
    print(f"{len(keys)} entries ({len(done)} already done), {a.workers} workers, shards from {shard_base}", flush=True)
    man_p = out / "manifest.tsv"
    cols = None
    t0, n_done, n_units, n_fail = time.time(), 0, 0, 0
    counter = mp.Value("i", 0)
    with mp.get_context("fork").Pool(a.workers, init_worker, (counter, out, a.antibody_cache, shard_base)) as pool, \
         open(done_p, "a") as fd, open(out / "failed.tsv", "a") as ff:
        mf = open(man_p, "a")
        for rows, dn, fl in pool.imap_unordered(run_chunk, chunks):
            if rows:
                if cols is None:
                    cols = list(rows[0])
                    if man_p.stat().st_size == 0:
                        mf.write("\t".join(cols) + "\n")
                for r in rows:
                    mf.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
            for k, n, sk in dn:
                fd.write(f"{k}\t{n}\t{sk}\n")
            for k, msg in fl:
                ff.write(f"{k}\t{msg}\n")
            mf.flush(); fd.flush(); ff.flush()
            n_done += len(dn); n_units += len(rows); n_fail += len(fl)
            if (n_done + n_fail) % 400 < a.chunk:
                el = time.time() - t0
                print(f"{n_done + n_fail}/{len(keys)} entries, {n_units} units, {n_fail} failed, {el:.0f}s "
                      f"({(n_done + n_fail) / max(el, 1):.1f}/s)", flush=True)
        mf.close()
    print(f"finished: {n_done} entries ok, {n_fail} failed, {n_units} units, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
