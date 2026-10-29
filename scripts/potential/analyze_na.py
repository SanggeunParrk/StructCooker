"""Protein-NA and NA-NA distance distributions from hist.npz -> JSON for the viewer.

usage: analyze_na.py <hist.npz> <out.json>
Protein-NA: 20 amino acids x 8 bases (DA DC DG DT RA RC RG RU), between chains. NA-NA: the 8 bases, 36 unordered
pairs, inside a chain and between chains. Unknown types (X, DN, RN) are left out. Per set: normalised distribution
and stats of every cell (raw or weighted counts), the pooled marginal, distributions pooled by amino acid and by
base, and Jensen-Shannon divergence (bits) to the marginal, with its shape-only part and the sampling-noise level.
"""
import json
import sys

import numpy as np

z = np.load(sys.argv[1])
T = [str(t) for t in z["types"]]
AA = T[:20]
NA = ["DA", "DC", "DG", "DT", "RA", "RC", "RG", "RU"]
NB = z["hist"].shape[-1] - 1
centres = (np.arange(NB) + 0.5) * float(z["bin"])
rng = np.random.default_rng(0)
ti = {t: T.index(t) for t in T}


def jsd(p, q):
    m = 0.5 * (p + q)
    kl = lambda a, b: float((a[a > 0] * np.log2(a[a > 0] / b[a > 0])).sum())  # noqa: E731
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def mom(q):
    mu = (q * centres).sum()
    return mu, np.sqrt((q * (centres - mu) ** 2).sum())


def shape_jsd(p, marg):
    mu, sd = mom(p)
    mum, sdm = mom(marg)
    if sd == 0 or sdm == 0:
        return jsd(p, marg)
    q = np.interp(mum + (centres - mu) * sdm / sd, centres, marg, left=0, right=0)
    return jsd(p, q / q.sum())


def near_jsd(p, marg, cut=8.0):
    k = centres < cut
    if p[k].sum() <= 0 or marg[k].sum() <= 0:     # a cell with nothing closer than cut
        return 0.0
    return jsd(p[k] / p[k].sum(), marg[k] / marg[k].sum())


def cell(H, s, a, b):
    i, j = sorted((ti[a], ti[b]))
    v = H[s, i, j, :NB].copy()
    return v + (H[s, j, i, :NB] if i != j else 0)


def stats_of(counts, marg, raw_n):
    n = counts.sum()
    if n <= 0:
        return None
    p = counts / n
    mu, sd = mom(p)
    null = np.mean([jsd(rng.multinomial(int(raw_n), marg) / raw_n, marg) for _ in range(20)]) if raw_n > 0 else 0
    return {"n": float(raw_n), "mean": round(float(mu), 3), "sd": round(float(sd), 3), "mode": float(centres[int(np.argmax(p))]),
            "jsd_marg": round(jsd(p, marg), 5), "jsd_shape": round(shape_jsd(p, marg), 5), "jsd_near": round(near_jsd(p, marg), 5), "null": round(float(null), 6)}


def block(rows, cols, s, sym):
    """Distributions of every (row, col) cell for scope s: raw and weighted."""
    out = {}
    for mode, key in (("raw", "hist"), ("weighted", "hist_weighted")):
        H, R = z[key], z["hist"]
        keys = [(a, b) for a in rows for b in cols if (not sym or rows.index(a) <= cols.index(b))]
        cells = {k: cell(H, s, *k) for k in keys}
        raw = {k: cell(R, s, *k) for k in keys}
        tot = sum(cells.values())
        marg = tot / tot.sum()
        d = {"marginal": [round(float(x), 5) for x in marg], "cells": [], "dist": [], "stats": [], "w": []}
        for k in keys:
            st = stats_of(cells[k], marg, raw[k].sum())
            d["cells"].append(list(k))
            d["dist"].append([round(float(x), 5) for x in (cells[k] / cells[k].sum())] if cells[k].sum() > 0 else [0] * NB)
            d["stats"].append(st)
            d["w"].append(float(cells[k].sum()))
        if not sym:   # protein x base: pooled by row and by column, and which factor moves the distribution more
            by_r = {a: sum(cells[(a, b)] for b in cols) for a in rows}
            by_c = {b: sum(cells[(a, b)] for a in rows) for b in cols}
            d["by_row"] = [[round(float(x), 5) for x in by_r[a] / by_r[a].sum()] for a in rows]
            d["by_col"] = [[round(float(x), 5) for x in by_c[b] / by_c[b].sum()] for b in cols]
            wr = np.array([by_r[a].sum() for a in rows]); wc = np.array([by_c[b].sum() for b in cols])
            jr = [jsd(by_r[a] / by_r[a].sum(), marg) for a in rows]; jc = [jsd(by_c[b] / by_c[b].sum(), marg) for b in cols]
            d["effect_row"] = round(float(np.average(jr, weights=wr)), 5)
            d["effect_col"] = round(float(np.average(jc, weights=wc)), 5)
            dna = sum(by_c[b] for b in cols[:4]); rna = sum(by_c[b] for b in cols[4:])
            d["dna_vs_rna"] = round(jsd(dna / dna.sum(), rna / rna.sum()), 5)
            d["mean_dna"], d["mean_rna"] = round(float((dna / dna.sum() * centres).sum()), 3), round(float((rna / rna.sum() * centres).sum()), 3)
        out[mode] = d
    return out


res = {"centres": centres.tolist(), "aa": AA, "na": NA,
       "pna": block(AA, NA, 1, False),
       "nana": {"intra": block(NA, NA, 0, True), "inter": block(NA, NA, 1, True)}}
json.dump(res, open(sys.argv[2], "w"), separators=(",", ":"))
print("written", sys.argv[2])
