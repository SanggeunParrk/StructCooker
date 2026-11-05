"""Reduce frames.npz to a small payload for the 3D viewer (sbatch only; needs ~1 GB of memory).

usage: frames_reduce.py <frames.npz> <skeleton.json> <out.json>
Works for both definitions: all-heavy-atom contacts (protein and bases) and side-chain contacts (protein only).
Dense grid (scope, reference type, partner type, 28^3 voxels of 1 A) -> 2 A voxels (14^3), partner types folded into 10
groups (all, AVILM, FWY, C, GP, STNQ, DE, KRH, DNA, RNA), reference types limited to the 20 amino acids and 8 bases.
Counts are stored as one byte per voxel on a log scale, decoded as expm1(q / 255 * log1p(max)) with the max kept per
(scope, reference, group). Both raw and weighted counts are written.
"""

import base64
import json
import sys

import numpy as np

z = np.load(sys.argv[1])
shape = tuple(int(x) for x in z["shape"])
types = [str(t) for t in z["types"]]
AA = types[:20]
SIDECHAIN = bool(z["sidechain"]) if "sidechain" in z.files else False
NA = (
    [] if SIDECHAIN else ["DA", "DC", "DG", "DT", "RA", "RC", "RG", "RU"]
)  # the side-chain definition is protein only
REFS = AA + NA
CLASSES = [
    ("AVILM", "AVILM"),
    ("FWY", "FWY"),
    ("CYS", "C"),
    ("GP", "GP"),
    ("STNQ", "STNQ"),
    ("DE", "DE"),
    ("KRH", "KRH"),
]
if (
    SIDECHAIN
):  # store every partner amino acid; the classes are sums of them, made in the page
    GROUPS = [("all", None)] + [(a, a) for a in AA]
    DISPLAY = ["all"] + [c for c, _ in CLASSES] + AA
    DERIVE = {c: list(m) for c, m in CLASSES}
    COLOR_GROUPS = AA
else:
    GROUPS = [
        ("all", None),
        *CLASSES,
        ("DNA", ["DA", "DC", "DG", "DT"]),
        ("RNA", ["RA", "RC", "RG", "RU"]),
    ]
    DISPLAY = [g for g, _ in GROUPS]
    DERIVE = {}
    COLOR_GROUPS = [c for c, _ in CLASSES] + ["DNA", "RNA"]
DEFINITION = (
    "곁사슬 heavy atom(N, CA, C, O 제외)끼리 5 Å 이내일 때만 접촉, 표준 아미노산만이며 Gly는 곁사슬이 없어 접촉이 없습니다"
    if SIDECHAIN
    else "heavy atom이 5 Å 이내이면 접촉, 아미노산과 핵산 모두"
)
members = {}
for name, spec in GROUPS:
    members[name] = [types.index(t) for t in (AA + NA if spec is None else list(spec))]
ref_idx = [types.index(t) for t in REFS]
NV = shape[-1] // 2


def dense(idx, val):
    d = np.zeros(int(np.prod(shape)), dtype=np.float32)
    d[idx] = val
    return d.reshape(shape)


out = {
    "key": "sc" if SIDECHAIN else "all",
    "definition": DEFINITION,
    "refs": REFS,
    "groups": DISPLAY,
    "stored": [g for g, _ in GROUPS],
    "derive": DERIVE,
    "colorGroups": COLOR_GROUPS,
    "scopes": [str(s) for s in z["scopes"]],
    "nv": NV,
    "vox": 2 * float(z["vox"]),
    "box": float(z["box"]),
    "skeleton": json.load(open(sys.argv[2])),
}
for mode, (ik, vk) in (("raw", ("idx", "val")), ("weighted", ("widx", "wval"))):
    D = dense(z[ik], z[vk])
    q_all, cmax, totals = [], [], []
    for s in range(2):
        for r in ref_idx:
            for g, _ in GROUPS:
                v = D[s, r][members[g]].sum(axis=0)  # (28, 28, 28)
                v = v.reshape(NV, 2, NV, 2, NV, 2).sum(axis=(1, 3, 5))  # 2 A voxels
                m = float(v.max())
                cmax.append(round(m, 3))
                totals.append(round(float(v.sum()), 1))
                q = (
                    np.zeros(v.shape, dtype=np.uint8)
                    if m <= 0
                    else np.round(255 * np.log1p(v) / np.log1p(m)).astype(np.uint8)
                )
                q_all.append(q.ravel())
    out[mode] = {
        "q": base64.b64encode(np.concatenate(q_all).tobytes()).decode(),
        "max": cmax,
        "total": totals,
    }
    print(mode, "payload", len(out[mode]["q"]) / 1e6, "MB base64")
json.dump(out, open(sys.argv[3], "w"), separators=(",", ":"))
print("written", sys.argv[3])
