"""Count residue-pair representative-atom distances over the training set (sbatch only).

usage: count_pairs.py --out DIR [--mode contact|all] [--rmax 15] [--workers 64] [--limit N] [--max-mb 30]
                      [--min-mb 0] [--min-sep 3] [--seed 0] [--weights entry_weights.tsv]
--mode contact: definition A (5 A heavy-atom contact), --mode all: definition B (every pair within --rmax).
Source: $BIOMOL_ROOT/lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb (one record per
PDB entry, waters removed). Per entry only the first assembly (lowest id), first model, altloc '.' (else first).
Writes DIR/hist.npz (hist and hist_weighted [scope, a, b, bin]; scopes intra / inter; weighted = each entry's counts
times its weight from --weights) and DIR/done.tsv.
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
ENV = None
MIN_SEP = 3
MODE = "contact"
RMAX_ALL = 15.0
WEIGHTS: dict[str, float] = {}


def init(min_sep, mode, rmax):
    global ENV, MIN_SEP, MODE, RMAX_ALL
    MIN_SEP, MODE, RMAX_ALL = min_sep, mode, rmax
    ENV = lmdb.open(
        str(DB), readonly=True, lock=False, readahead=False, max_readers=512
    )


def pick(keys):
    """Lowest assembly id, then lowest model, altloc '.' preferred."""

    def k(c):
        a, m, alt = c.split("_")
        return (int(a), int(m), alt != ".", alt)

    return min(keys, key=k)


def one(key):
    with ENV.begin() as txn:
        raw = txn.get(key.encode())
    rec = from_bytes(bytes(raw))
    ck = pick(rec)
    cm = CIFMolAttached.from_dict(rec[ck]["cifmol_attached_dict"])
    if MODE == "all":
        return pairs.count_assembly_all(cm, rmax=RMAX_ALL, min_sep=MIN_SEP)
    return pairs.count_assembly(cm, min_sep=MIN_SEP)


def chunk(keys):
    nb = int(RMAX_ALL / pairs.BIN) if MODE == "all" else pairs.NBIN
    total = np.zeros((2, pairs.T, pairs.T, nb + 1))
    total_w = np.zeros_like(total)
    done, fail = [], []
    for k in keys:
        try:
            h = one(k)
            total += h
            total_w += h * WEIGHTS.get(k, 1.0)
            done.append((k, int(h[0].sum()), int(h[1].sum())))
        except Exception:  # noqa: BLE001
            fail.append((k, traceback.format_exc().strip().splitlines()[-1][:200]))
    return total, total_w, done, fail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=64)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-mb", type=float, default=0)
    ap.add_argument("--max-mb", type=float, default=30)
    ap.add_argument("--min-sep", type=int, default=3)
    ap.add_argument(
        "--mode",
        choices=["contact", "all"],
        default="contact",
        help="contact = definition A (heavy atoms within 5 A), all = definition B (every pair within --rmax)",
    )
    ap.add_argument(
        "--rmax",
        type=float,
        default=15.0,
        help="definition B: distance limit between representative atoms",
    )
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=8)
    ap.add_argument(
        "--weights",
        default=None,
        help="entry_weights.tsv (entry, n, weight); default weight 1",
    )
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for line in open(str(DB) + ".index.tsv"):
        if line.startswith(("#", "key")):
            continue
        k, b = line.split("\t")[:2]
        sizes[k] = int(b)
    if a.weights:
        WEIGHTS.update(
            {l.split("\t")[0]: float(l.split("\t")[2]) for l in open(a.weights)}
        )
    keys = [k for k, b in sizes.items() if a.min_mb * 1e6 <= b < a.max_mb * 1e6]
    random.Random(a.seed).shuffle(keys)
    if a.limit:
        keys = keys[: a.limit]
    keys.sort(key=lambda k: -sizes[k])
    chunks = [keys[i : i + a.chunk] for i in range(0, len(keys), a.chunk)]
    print(f"{len(keys)} entries, {a.workers} workers", flush=True)
    nbins = int(a.rmax / pairs.BIN) if a.mode == "all" else pairs.NBIN
    hist = np.zeros((2, pairs.T, pairs.T, nbins + 1))
    hist_w = np.zeros_like(hist)
    n_ok = n_fail = 0
    t0 = time.time()
    with (
        mp.get_context("fork").Pool(
            a.workers, init, (a.min_sep, a.mode, a.rmax)
        ) as pool,
        open(out / "done.tsv", "w") as fd,
        open(out / "failed.tsv", "w") as ff,
    ):
        for h, hw, dn, fl in pool.imap_unordered(chunk, chunks):
            hist += h
            hist_w += hw
            n_ok += len(dn)
            n_fail += len(fl)
            for k, ni, nx in dn:
                fd.write(f"{k}\t{ni}\t{nx}\n")
            for k, m in fl:
                ff.write(f"{k}\t{m}\n")
            if (n_ok + n_fail) % 400 < a.chunk:
                print(
                    f"{n_ok + n_fail}/{len(keys)} ({n_fail} failed) {time.time() - t0:.0f}s",
                    flush=True,
                )
    np.savez_compressed(
        out / "hist.npz",
        hist=hist,
        hist_weighted=hist_w,
        types=np.array(pairs.TYPES),
        bin=pairs.BIN,
        rmax=a.rmax if a.mode == "all" else pairs.RMAX,
        mode=a.mode,
        scopes=np.array(pairs.SCOPES),
        min_sep=a.min_sep,
    )
    print(
        f"done: {n_ok} entries, {n_fail} failed; pairs intra {hist[0].sum():.0f} inter {hist[1].sum():.0f}; "
        f"{time.time() - t0:.0f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
