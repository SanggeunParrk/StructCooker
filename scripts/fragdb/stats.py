"""Sanity checks and size/redundancy statistics of fragment-pair parts (sbatch only). usage: stats.py <dir> [k=7]"""

import datetime
import sys
from pathlib import Path

import numpy as np
import zstandard

d = Path(sys.argv[1])
k = int(sys.argv[2]) if len(sys.argv) > 2 else 7
parts = sorted(p for p in d.glob("part_*") if p.is_dir())
fields = (
    "seq",
    "coords",
    "center",
    "chain",
    "cluster",
    "inter",
    "n_atom_pairs",
    "pdb",
    "assembly",
    "date",
    "resolution",
)
A = {f: np.concatenate([np.load(p / f"{f}.npy") for p in parts]) for f in fields}
A["pdb"] = A["pdb"].astype("<u4").view("S4").astype("U4")  # unpack the 4-character ids
n = len(A["inter"])
c = A["coords"].astype(np.float32) / 100.0  # (n, 2, k, 5, 3)
print(
    f"{n:,} pairs from {len(set(A['pdb'].tolist())):,} entries; inter-chain {A['inter'].mean():.0%}, intra {1 - A['inter'].mean():.0%}"
)
print(
    f"bytes per pair as stored: {sum(a.nbytes for a in A.values()) / n:.0f}  (coordinates alone {A['coords'].nbytes / n:.0f})"
)
mid = k // 2
print(
    "\ncentre residue of fragment A in its own frame (must be the same everywhere): mean +- sd"
)
for i, a in enumerate(("N", "CA", "C", "O", "CB")):
    v = c[:, 0, mid, i]
    print(
        f"  {a:2s}",
        "  ".join(
            f"{m:6.2f}+-{s:.2f}" for m, s in zip(v.mean(0), v.std(0), strict=False)
        ),
    )
pb = np.linalg.norm(c[:, :, :-1, 2] - c[:, :, 1:, 0], axis=-1)  # C(i) - N(i+1)
print(
    f"\npeptide bond C-N inside the fragments: mean {pb.mean():.3f} A, max {pb.max():.2f} A (expect ~1.33)"
)
dcb = np.linalg.norm(c[:, 0, mid, 4] - c[:, 1, mid, 4], axis=-1)
print(
    f"Cb-Cb distance between the two centres: mean {dcb.mean():.2f} A, 5-95%: {np.percentile(dcb, 5):.1f}-{np.percentile(dcb, 95):.1f} A, max {dcb.max():.1f} A"
)
dt = A["date"]
f = lambda x: str(datetime.date(1970, 1, 1) + datetime.timedelta(days=int(x)))
print(
    f"release dates {f(dt.min())} .. {f(dt.max())}; pairs without a date: {(dt == 0).sum()}; resolution median {np.nanmedian(A['resolution'].astype(np.float32)):.2f} A"
)
# exact sequence repeats
seq = A["seq"].astype(np.uint64)


def pack(s):  # 7 letters of 5 bits
    v = np.zeros(len(s), dtype=np.uint64)
    for i in range(s.shape[1]):
        v = (v << np.uint64(5)) | s[:, i]
    return v


ka, kb = pack(seq[:, 0]), pack(seq[:, 1])
lo, hi = np.minimum(ka, kb), np.maximum(ka, kb)
pair = lo * np.uint64(1 << 20) ^ hi * np.uint64(
    0x9E3779B97F4A7C15
)  # hash of the unordered sequence pair
u, cnt = np.unique(pair, return_counts=True)
print(
    f"\nexact sequence pairs: {len(u):,} distinct among {n:,} ({100 * len(u) / n:.1f}%); seen more than once: {(cnt > 1).sum():,} types covering {cnt[cnt > 1].sum() / n:.1%} of the pairs"
)
ua = np.unique(ka)
print(
    f"distinct 7-mers (fragment A): {len(ua):,}; share of the 20^7 = 1.28e9 possible: {len(ua) / 20**7:.2e}"
)
# near-duplicates: pairs of the same entry and chain pair with centres within 2 residues on both sides
key = np.char.add(np.char.add(A["pdb"], A["chain"][:, 0]), A["chain"][:, 1])
order = np.lexsort((A["center"][:, 1], A["center"][:, 0], key))
kk, ca, cb = key[order], A["center"][order, 0], A["center"][order, 1]
near = np.zeros(n, dtype=bool)
for s in range(1, 6):
    same = kk[s:] == kk[:-s]
    close = same & (np.abs(ca[s:] - ca[:-s]) <= 2) & (np.abs(cb[s:] - cb[:-s]) <= 2)
    near[s:] |= close
print(
    f"\nredundancy: {near.mean():.0%} of the pairs have a neighbouring pair (centres within 2 residues on both sides) in the same chain pair"
)
# clusters
cl = A["cluster"]
same_cl = (cl[:, 0] == cl[:, 1]) & (cl[:, 0] > 0)
print(f"pairs whose two fragments come from the same 30% cluster: {same_cl.mean():.0%}")
raw = np.ascontiguousarray(A["coords"][: min(n, 200_000)]).tobytes()
print(
    f"zstd-3 of the coordinates: {len(zstandard.ZstdCompressor(level=3).compress(raw)) / len(raw):.2f} of the raw size"
)
