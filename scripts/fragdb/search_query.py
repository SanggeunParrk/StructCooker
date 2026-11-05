"""Search the fragment-pair DB with the 7-mer pairs of one query structure (sbatch only).

usage: search_query.py --pdb 9yc5 --roles A_1:antigen,B_1:H,C_1:L --db DIR --out DIR --k 7 [--assembly 1] [--topk N] [--same-kind]
The query's pairs are extracted exactly as the DB was built (structcooker.instructions.fragdb.extract). For every query pair the DB entries are scored
by sequence identity over the 14 positions (both orders of the two fragments are tried), the best ``topk`` are compared with the query
by frame RMSD: both pairs are expressed in the frame of the first fragment's centre residue, so no superposition is needed. RMSD is over
all stored atoms (N CA C O CB of 14 residues) and, separately, over the second fragment only (its placement relative to the first).
A random sample of DB pairs gives the RMSD one would get without any sequence information. With ``--same-kind`` a query
between two chains is compared only with DB pairs between two chains, a query inside one chain only with DB pairs inside one chain.
"""

import argparse
import datetime
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import lmdb
import numpy as np

from structcooker.instructions.fragdb.extract import extract_pairs
from structcooker.instructions.mols import CIFMolAttached
from structcooker.instructions.potential.pairs import AA
from structcooker.instructions.transforms.codecs import from_bytes

ROOT = Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean"))
K, MID, MINID = 7, 3, 8
S = {}


def load(db):
    parts = sorted(p for p in Path(db).glob("part_*") if p.is_dir())
    S["parts"] = parts
    S["seq"] = [np.load(p / "seq.npy") for p in parts]
    sizes = [len(x) for x in S["seq"]]
    S["off"] = np.concatenate([[0], np.cumsum(sizes)])
    seq = np.concatenate(S["seq"])
    S["N"] = len(seq)
    S["A"] = np.ascontiguousarray(
        seq[:, 0].T
    )  # (7, N): column-wise scans are cache friendly
    S["B"] = np.ascontiguousarray(seq[:, 1].T)
    del seq
    for f in ("pdb", "date", "center", "inter", "chain", "resolution"):
        S[f] = np.concatenate([np.load(p / f"{f}.npy") for p in parts])
    S["coords"] = [np.load(p / "coords.npy", mmap_mode="r") for p in parts]


def coords_of(idx):
    part = np.searchsorted(S["off"], idx, side="right") - 1
    out = np.empty((len(idx), 2, K, 5, 3), dtype=np.float32)
    for pi in np.unique(part):
        m = part == pi
        out[m] = S["coords"][pi][idx[m] - S["off"][pi]].astype(np.float32) / 100.0
    return out


def to_b_frame(c):
    """Re-express stored pairs (m, 2, k, 5, 3), given in A's frame, in the frame of B's centre residue, B first."""
    n, ca, cc = c[:, 1, MID, 0], c[:, 1, MID, 1], c[:, 1, MID, 2]
    e1 = cc - ca
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    v = n - ca
    e2 = v - (v * e1).sum(1, keepdims=True) * e1
    e2 /= np.linalg.norm(e2, axis=1, keepdims=True)
    r = np.stack([e1, e2, np.cross(e1, e2)], axis=2)
    loc = np.einsum("mji,mfkaj->mfkai", r, c - ca[:, None, None, None, :])
    return loc[:, ::-1]


def identity(qa, qb):
    s1 = np.zeros(S["N"], dtype=np.uint8)
    s2 = np.zeros(S["N"], dtype=np.uint8)
    for i in range(K):
        s1 += S["A"][i] == qa[i]
        s1 += S["B"][i] == qb[i]
        s2 += S["A"][i] == qb[i]
        s2 += S["B"][i] == qa[i]
    return s1, s2


def one_query(args):
    qi, qseq, qcoord, topk, want = args
    np.random.default_rng(qi)
    s1, s2 = identity(qseq[0], qseq[1])
    score = np.maximum(s1, s2)
    if want >= 0:  # same kind of pair only: between chains, or inside one chain
        score = np.where(S["inter"] == want, score, 0).astype(np.uint8)
    top = np.argpartition(score, -topk)[-topk:]
    top = top[np.argsort(-score[top], kind="stable")]
    swapped = s2[top] > s1[top]
    c = coords_of(top)
    if swapped.any():
        c[swapped] = to_b_frame(c[swapped])
    d = np.linalg.norm(c - qcoord[None], axis=-1)  # (topk, 2, k, 5)
    rmsd_all = np.sqrt((d**2).mean(axis=(1, 2, 3)))
    rmsd_b = np.sqrt((d[:, 1] ** 2).mean(axis=(1, 2)))
    return qi, int(score.max()), top, score[top], swapped, rmsd_all, rmsd_b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument(
        "--roles",
        required=True,
        help="asym_oper chain ids and roles, e.g. A_1:antigen,B_1:H,C_1:L",
    )
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--assembly", default="1")
    ap.add_argument("--k", type=int, default=7)
    ap.add_argument(
        "--min-id",
        type=int,
        default=0,
        help="identity (of 2k positions) that counts as a sequence match; default ~57%%",
    )
    ap.add_argument(
        "--topk",
        type=int,
        default=0,
        help="candidates compared geometrically per query pair (default 200 / 500 / 2000 for k = 7 / 5 / 3)",
    )
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--max-per-class", type=int, default=400)
    ap.add_argument(
        "--same-kind",
        action="store_true",
        help="a between-chains query is matched only with between-chains DB pairs (and likewise inside one chain)",
    )
    a = ap.parse_args()
    global K, MID, MINID
    K, MID = a.k, a.k // 2
    MINID = a.min_id or int(np.ceil(0.57 * 2 * K))
    a.topk = a.topk or {7: 200, 5: 500, 3: 2000}.get(K, 200)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    roles = dict(x.split(":") for x in a.roles.split(","))
    t0 = time.time()
    env = lmdb.open(
        str(ROOT / "lmdb/pdb/cif_attached/cif_pdb_attached.lmdb"),
        readonly=True,
        lock=False,
    )
    with env.begin() as txn:
        rec = from_bytes(bytes(txn.get(a.pdb.encode())))
    key = min(
        (k for k in rec if k.split("_")[0] == a.assembly),
        key=lambda c: (int(c.split("_")[1]), c.split("_")[2] != ".", c),
    )
    cm = CIFMolAttached.from_dict(rec[key]["cifmol_attached_dict"])
    res = extract_pairs(cm, k=K)
    chains = res["chain"]
    keep = np.array([c[0] in roles and c[1] in roles for c in chains])
    cls = []
    for i in np.flatnonzero(keep):
        r1, r2 = roles[chains[i][0]], roles[chains[i][1]]
        inter = bool(res["inter"][i])
        if inter and "antigen" in (r1, r2):
            cls.append("Ab-Ag interface")
        elif inter:
            cls.append("H-L interface")
        else:
            cls.append("inside " + ("antigen" if r1 == "antigen" else "antibody"))
    idx = np.flatnonzero(keep)
    cls = np.array(cls)
    rng = np.random.default_rng(0)
    sel = []
    for c in sorted(set(cls)):
        ii = idx[cls == c]
        sel += (
            list(rng.choice(ii, min(len(ii), a.max_per_class), replace=False))
            if len(ii) > a.max_per_class
            else list(ii)
        )
    sel = np.array(sorted(sel))
    cls_of = dict(zip(idx.tolist(), cls.tolist(), strict=False))
    print(
        f"query {a.pdb} assembly {a.assembly}, released: see metadata; {len(res['inter'])} pairs in the structure, {keep.sum()} among the three chains "
        f"{ {c: int((cls == c).sum()) for c in sorted(set(cls))} }; searching {len(sel)}",
        flush=True,
    )
    load(a.db)
    print(f"DB loaded: {S['N']:,} pairs in {time.time() - t0:.0f}s", flush=True)
    qseq = res["seq"][sel]
    qcoord = res["coords"][sel].astype(np.float32) / 100.0

    # null: RMSD of random DB pairs against each query pair
    def null_for(
        kind,
    ):  # RMSD of random DB pairs (of one kind) against every query pair
        pool_idx = (
            np.flatnonzero(S["inter"] == kind) if a.same_kind else np.arange(S["N"])
        )
        rc = coords_of(np.sort(rng.choice(pool_idx, 4000, replace=False)))
        return np.sqrt(
            ((rc[None, :, 1] - qcoord[:, None, 1]) ** 2).sum(-1).mean(axis=(2, 3))
        )

    null_by = (
        {0: null_for(0), 1: null_for(1)}
        if a.same_kind
        else {0: null_for(0), 1: null_for(0)}
    )
    qkind = res["inter"][sel].astype(int)
    null_b = np.where(qkind[:, None] == 1, null_by[1], null_by[0])
    with mp.get_context("fork").Pool(a.workers) as pool:
        results = pool.map(
            one_query,
            [
                (
                    i,
                    qseq[i],
                    qcoord[i],
                    a.topk,
                    int(res["inter"][sel[i]]) if a.same_kind else -1,
                )
                for i in range(len(sel))
            ],
        )
    print(f"searched in {time.time() - t0:.0f}s", flush=True)
    rows = []
    rng_null = np.random.default_rng(1)
    for qi, best, top, sc, sw, ra, rb in results:
        j = int(np.argmin(rb))
        i0 = 0
        rows.append(
            {
                "query": int(sel[qi]),
                "class": cls_of[int(sel[qi])],
                "chains": list(map(str, chains[sel[qi]])),
                "center": list(map(int, res["center"][sel[qi]])),
                "seqA": "".join(AA[x] for x in qseq[qi, 0]),
                "seqB": "".join(AA[x] for x in qseq[qi, 1]),
                "best_id": best,
                "top1_rmsd_all": float(ra[i0]),
                "top1_rmsd_b": float(rb[i0]),
                "min_rmsd_b": float(rb[j]),
                "min_rmsd_b_id": int(sc[j]),
                "n_hits_id8": int((sc >= MINID).sum()),
                "n_good": int(((sc >= MINID) & (rb <= 2.5)).sum()),
                "null_min_rmsd_b": float(
                    np.median(
                        [
                            np.min(rng_null.choice(null_b[qi], a.topk))
                            for _ in range(100)
                        ]
                    )
                ),
                "null_rmsd_b_median": float(np.median(null_b[qi])),
                "null_frac_le2.5": float((null_b[qi] <= 2.5).mean()),
                "top1": {
                    "pdb": int(S["pdb"][top[i0]]).to_bytes(4, "little").decode(),
                    "chains": list(map(str, S["chain"][top[i0]])),
                    "center": list(map(int, S["center"][top[i0]])),
                    "date": str(
                        datetime.date(1970, 1, 1)
                        + datetime.timedelta(days=int(S["date"][top[i0]]))
                    ),
                    "swapped": bool(sw[i0]),
                },
            }
        )
    (out / "results.json").write_text(
        json.dumps(
            {
                "pdb": a.pdb,
                "k": K,
                "min_id": MINID,
                "topk": a.topk,
                "roles": roles,
                "rows": rows,
                "null_b": np.percentile(null_b, [1, 5, 25, 50, 75]).tolist(),
                "db_pairs": S["N"],
            }
        )
    )
    print("written", out / "results.json")


if __name__ == "__main__":
    main()
