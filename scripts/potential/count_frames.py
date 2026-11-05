"""Contacts in every residue's backbone frame, pooled by residue type (sbatch only).

usage: count_frames.py --out DIR [--workers 56] [--limit N] [--max-mb 1e6] [--seed 0] [--weights entry_weights.tsv]
                       [--sidechain]
Definition A contacts (heavy atoms within 5 A) of the training set; the partner's representative atom is placed in the
frame of the other residue. A grid of 1 A voxels over +-14 A, indexed [scope, reference type, partner type, x, y, z]
(see structcooker.instructions.potential.pairs.count_frame_contacts). Each worker takes a fixed share of the entries and writes
DIR/part_NNN.npz; the parts are summed into DIR/frames.npz (non-zero cells only: idx, val, widx, wval).
"""

import argparse
import os

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"

import multiprocessing as mp  # noqa: E402
import random  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import lmdb  # noqa: E402
import numpy as np  # noqa: E402

from structcooker.instructions.mols import CIFMolAttached  # noqa: E402
from structcooker.instructions.potential import pairs  # noqa: E402
from structcooker.instructions.transforms.codecs import from_bytes  # noqa: E402

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
DB = ROOT / "lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb"
SHAPE = (2, pairs.T, pairs.T, pairs.NV, pairs.NV, pairs.NV)
SIDECHAIN = False
WEIGHTS: dict[str, float] = {}


def pick(keys):
    def k(c):
        a, m, alt = c.split("_")
        return (int(a), int(m), alt != ".", alt)

    return min(keys, key=k)


def run_shard(task):
    sid, keys, out = task
    env = lmdb.open(
        str(DB), readonly=True, lock=False, readahead=False, max_readers=512
    )
    grid = np.zeros(SHAPE, dtype=np.float32)
    wgrid = np.zeros(SHAPE, dtype=np.float32)
    n_in = n_out = n_fail = 0
    fails = []
    t0 = time.time()
    for k, key in enumerate(keys):
        try:
            with env.begin() as txn:
                raw = txn.get(key.encode())
            rec = from_bytes(bytes(raw))
            cm = CIFMolAttached.from_dict(rec[pick(rec)]["cifmol_attached_dict"])
            _, a, b = pairs.count_frame_contacts(
                cm,
                grid=grid,
                wgrid=wgrid,
                weight=WEIGHTS.get(key, 1.0),
                sidechain=SIDECHAIN,
            )
            n_in += a
            n_out += b
        except Exception:  # noqa: BLE001
            n_fail += 1
            fails.append((key, traceback.format_exc().strip().splitlines()[-1][:160]))
        if (k + 1) % 1000 == 0:
            print(
                f"shard {sid}: {k + 1}/{len(keys)} entries {time.time() - t0:.0f}s",
                flush=True,
            )
    nz = np.flatnonzero(grid.ravel())
    nzw = np.flatnonzero(wgrid.ravel())
    np.savez_compressed(
        Path(out) / f"part_{sid:03d}.npz",
        idx=nz,
        val=grid.ravel()[nz],
        widx=nzw,
        wval=wgrid.ravel()[nzw],
        n_in=n_in,
        n_out=n_out,
    )
    return sid, len(keys), n_in, n_out, n_fail, fails[:3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=56)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--max-mb", type=float, default=1e6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--weights", default=None)
    ap.add_argument(
        "--sidechain",
        action="store_true",
        help="narrow contact: standard amino acids only, side-chain heavy atoms within 5 A of each other",
    )
    a = ap.parse_args()
    global SHAPE, SIDECHAIN
    SIDECHAIN = a.sidechain
    if a.sidechain:
        SHAPE = (2, len(pairs.AA), len(pairs.AA), pairs.NV, pairs.NV, pairs.NV)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.weights:
        WEIGHTS.update(
            {
                line.split("\t")[0]: float(line.split("\t")[2])
                for line in open(a.weights)
            }
        )
    sizes = {}
    for line in open(str(DB) + ".index.tsv"):
        if line.startswith(("#", "key")):
            continue
        k, b = line.split("\t")[:2]
        sizes[k] = int(b)
    keys = [k for k, b in sizes.items() if b < a.max_mb * 1e6]
    random.Random(a.seed).shuffle(keys)
    if a.limit:
        keys = keys[: a.limit]
    keys.sort(key=lambda k: -sizes[k])
    tasks = [(s, keys[s :: a.workers], str(out)) for s in range(a.workers)]
    print(
        f"{len(keys)} entries, {a.workers} shards, grid {SHAPE} = {np.prod(SHAPE) * 4 / 1e6:.0f} MB per worker (x2)",
        flush=True,
    )
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as pool:
        res = list(pool.imap_unordered(run_shard, tasks))
    n_in = sum(r[2] for r in res)
    n_out = sum(r[3] for r in res)
    n_fail = sum(r[4] for r in res)
    print(
        f"workers done in {time.time() - t0:.0f}s: {n_in} contacts placed, {n_out} outside the box, {n_fail} entries failed",
        flush=True,
    )
    for r in res:
        for f in r[5]:
            print("  fail:", f, flush=True)
    acc, wacc = (
        np.zeros(int(np.prod(SHAPE)), dtype=np.float32),
        np.zeros(int(np.prod(SHAPE)), dtype=np.float32),
    )
    for s in range(a.workers):
        z = np.load(out / f"part_{s:03d}.npz")
        acc[z["idx"]] += z["val"]
        wacc[z["widx"]] += z["wval"]
    nz, nzw = np.flatnonzero(acc), np.flatnonzero(wacc)
    np.savez_compressed(
        out / "frames.npz",
        idx=nz,
        val=acc[nz],
        widx=nzw,
        wval=wacc[nzw],
        shape=np.array(SHAPE),
        box=pairs.BOX,
        vox=pairs.VOX,
        types=np.array(pairs.TYPES[: SHAPE[1]]),
        scopes=np.array(pairs.SCOPES),
        n_in=n_in,
        n_out=n_out,
        sidechain=a.sidechain,
    )
    for s in range(a.workers):
        (out / f"part_{s:03d}.npz").unlink()
    print(
        f"frames.npz: {len(nz)} non-zero cells of {acc.size}, {os.path.getsize(out / 'frames.npz') / 1e6:.0f} MB; total {time.time() - t0:.0f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
