"""Statistics of a full fragment-pair DB, part by part (sbatch only). usage: stats_full.py <dir> [k=7]

Coordinate checks use the first 400k pairs of the first part; counts, dates, redundancy (per part: an entry never spans two parts)
and the exact-sequence repeat rate use all pairs.
"""

import datetime
import sys
from pathlib import Path

import numpy as np
import zstandard

d = Path(sys.argv[1])
k = int(sys.argv[2]) if len(sys.argv) > 2 else 7
parts = sorted(p for p in d.glob("part_*") if p.is_dir())
P = np.uint64(0x9E3779B97F4A7C15)
n_pairs = n_inter = n_near = 0
dmin, dmax, ndate0 = 10**9, 0, 0
nbytes = 0
res_all, entries = [], set()
hashes, ka_all = [], []
for p in parts:
    seq = np.load(p / "seq.npy").astype(np.uint64)
    n = len(seq)
    n_pairs += n
    n_inter += int(np.load(p / "inter.npy").sum())
    nbytes += sum(f.stat().st_size for f in p.glob("*.npy"))
    date = np.load(p / "date.npy")
    dmin, dmax, ndate0 = (
        min(dmin, int(date.min())),
        max(dmax, int(date.max())),
        ndate0 + int((date == 0).sum()),
    )
    res_all.append(
        np.load(p / "resolution.npy").astype(np.float32)[:: max(1, n // 20000)]
    )
    pdb = np.load(p / "pdb.npy")
    entries.update(np.unique(pdb).tolist())

    def pack(s):
        v = np.zeros(len(s), dtype=np.uint64)
        for i in range(s.shape[1]):
            v = (v << np.uint64(5)) | s[:, i]
        return v

    ka, kb = pack(seq[:, 0]), pack(seq[:, 1])
    lo, hi = np.minimum(ka, kb), np.maximum(ka, kb)
    hashes.append(lo * P + hi * np.uint64(0xC2B2AE3D27D4EB4F))
    ka_all.append(np.unique(ka))  # mod 2^64, one value per unordered sequence pair
    ch = np.load(p / "chain.npy").view(np.uint32).reshape(n, 2, -1)
    chash = (ch.astype(np.uint64) * (np.arange(ch.shape[-1], dtype=np.uint64) + 1)).sum(
        -1
    )
    cen = np.load(p / "center.npy")
    order = np.lexsort((cen[:, 1], cen[:, 0], chash[:, 1], chash[:, 0], pdb))
    kk = np.stack(
        [
            pdb[order].astype(np.int64),
            chash[order, 0].astype(np.int64),
            chash[order, 1].astype(np.int64),
        ],
        axis=1,
    )
    ca, cb = cen[order, 0], cen[order, 1]
    near = np.zeros(n, dtype=bool)
    for s in range(1, 6):
        same = (kk[s:] == kk[:-s]).all(axis=1)
        near[s:] |= (
            same & (np.abs(ca[s:] - ca[:-s]) <= 2) & (np.abs(cb[s:] - cb[:-s]) <= 2)
        )
    n_near += int(near.sum())
    print(f"{p.name}: {n:,} pairs", flush=True)
f = lambda x: str(datetime.date(1970, 1, 1) + datetime.timedelta(days=int(x)))
print(
    f"\n{n_pairs:,} pairs from {len(entries):,} entries in {len(parts)} parts; inter-chain {n_inter / n_pairs:.1%}"
)
print(
    f"on disk: {nbytes / 1e9:.1f} GB ({nbytes / n_pairs:.0f} B per pair); release dates {f(dmin)} .. {f(dmax)}; pairs without a date: {ndate0}"
)
print(f"resolution median {np.nanmedian(np.concatenate(res_all)):.2f} A")
h = np.concatenate(hashes)
u, cnt = np.unique(h, return_counts=True)
print(
    f"exact sequence pairs: {len(u):,} distinct of {n_pairs:,} ({100 * len(u) / n_pairs:.1f}%); types seen more than once cover {cnt[cnt > 1].sum() / n_pairs:.1%} of the pairs; the most frequent type occurs {cnt.max():,} times"
)
ua = np.unique(np.concatenate(ka_all))
print(f"distinct 7-mers (fragment A): {len(ua):,} ({len(ua) / 20**k:.1e} of 20^{k})")
print(
    f"redundancy: {n_near / n_pairs:.0%} of the pairs have a neighbouring pair (centres within 2 residues on both sides) in the same chain pair"
)
# coordinate checks on a sample
p0 = parts[0]
m = 400_000
c = np.load(p0 / "coords.npy", mmap_mode="r")[:m].astype(np.float32) / 100.0
mid = k // 2
print("\ncentre residue of A in its frame (mean +- sd):")
for i, a in enumerate(("N", "CA", "C", "O", "CB")):
    v = c[:, 0, mid, i]
    print(
        f"  {a:2s}",
        "  ".join(
            f"{x:6.2f}+-{s:.2f}" for x, s in zip(v.mean(0), v.std(0), strict=False)
        ),
    )
pb = np.linalg.norm(c[:, :, :-1, 2] - c[:, :, 1:, 0], axis=-1)
print(
    f"peptide bond C-N: mean {pb.mean():.3f} A, max {pb.max():.2f} A; centre Cb-Cb: mean {np.linalg.norm(c[:, 0, mid, 4] - c[:, 1, mid, 4], axis=-1).mean():.2f} A"
)
raw = np.ascontiguousarray(
    np.load(p0 / "coords.npy", mmap_mode="r")[:200_000]
).tobytes()
print(
    f"zstd-3 of the coordinates: {len(zstandard.ZstdCompressor(level=3).compress(raw)) / len(raw):.2f} of the raw size"
)
