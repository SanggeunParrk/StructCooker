"""Compare the residue-pair distance distributions of hist.npz (amino acids only), write JSON for the viewer.

usage: analyze_hist.py <hist.npz> <out.json>
Per scope (intra / inter) and weighting (raw / weighted): the normalised distribution of every unordered pair of the
20 standard amino acids, the all-pair marginal, Jensen-Shannon divergence (bits) of each pair to the marginal and
to its other-scope twin, the part of it that is shape rather than position and width (``jsd_shape``: the average pair
shifted and stretched to the pair's mean and sd), and the divergence a pair of that size would show by sampling noise alone (multinomial
draws from the marginal with the pair's raw count).
"""

import json
import sys

import numpy as np

z = np.load(sys.argv[1])
types = [str(t) for t in z["types"]]
AA = types[:20]
NB = (
    z["hist"].shape[-1] - 1
)  # 40 for definition A (0-20 A), 30 for definition B (0-15 A); the last bin is overflow
centres = (np.arange(NB) + 0.5) * float(z["bin"])
rng = np.random.default_rng(0)
iu = [(i, j) for i in range(20) for j in range(i, 20)]


def jsd(p, q):
    m = 0.5 * (p + q)

    def kl(a, b):
        ok = a > 0
        return float((a[ok] * np.log2(a[ok] / b[ok])).sum())

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def shape_jsd(p, marg):
    """Divergence left after the average pair is moved and stretched to the pair's mean and width."""

    def mom(q):
        mu = (q * centres).sum()
        return mu, np.sqrt((q * (centres - mu) ** 2).sum())

    mu, sd = mom(p)
    mum, sdm = mom(marg)
    q = np.interp(mum + (centres - mu) * sdm / sd, centres, marg, left=0, right=0)
    return jsd(p, q / q.sum())


def near_jsd(p, marg, cut=8.0):
    """Divergence of the two distributions restricted to r < cut and renormalised (the range where chemistry shows)."""
    k = centres < cut
    if p[k].sum() <= 0 or marg[k].sum() <= 0:  # a cell with nothing closer than cut
        return 0.0
    return jsd(p[k] / p[k].sum(), marg[k] / marg[k].sum())


def pair_hist(H, s):
    """{(i, j): counts over the 40 bins} for the unordered amino-acid pairs, overflow dropped."""
    return {(i, j): H[s, i, j, :NB] + (0 if i == j else H[s, j, i, :NB]) for i, j in iu}


out = {
    "aa": AA,
    "centres": centres.tolist(),
    "pairs": [f"{AA[i]}{AA[j]}" for i, j in iu],
    "scopes": {},
}
raw = {s: pair_hist(z["hist"], s) for s in (0, 1)}
wei = {s: pair_hist(z["hist_weighted"], s) for s in (0, 1)}
for s, sname in enumerate(("intra", "inter")):
    d = {}
    for mode, P in (("raw", raw[s]), ("weighted", wei[s])):
        tot = sum(P.values())
        marg = tot / tot.sum()
        dist, stats = [], []
        for k, (i, j) in enumerate(iu):
            h = P[(i, j)]
            p = h / h.sum()
            n_raw = raw[s][(i, j)].sum()
            null = np.mean(
                [
                    jsd(rng.multinomial(int(n_raw), marg) / n_raw, marg)
                    for _ in range(30)
                ]
            )
            mean = float((p * centres).sum())
            stats.append(
                {
                    "n": float(n_raw),
                    "mean": round(mean, 3),
                    "sd": round(float(np.sqrt((p * (centres - mean) ** 2).sum())), 3),
                    "mode": float(centres[int(np.argmax(p))]),
                    "jsd_marg": round(jsd(p, marg), 5),
                    "jsd_shape": round(shape_jsd(p, marg), 5),
                    "jsd_near": round(near_jsd(p, marg), 5),
                    "null": round(float(null), 6),
                    "share_lt5": round(float(p[:10].sum()), 4),
                    "share_gt9": round(float(p[18:].sum()), 4),
                }
            )
            dist.append([round(float(x), 5) for x in p])
        d[mode] = {
            "marginal": [round(float(x), 5) for x in marg],
            "dist": dist,
            "stats": stats,
        }
    out["scopes"][sname] = d
# intra vs inter
out["cross"] = {}
for mode, W in (("raw", raw), ("weighted", wei)):
    ti, tx = sum(W[0].values()), sum(W[1].values())
    mi, mx = ti / ti.sum(), tx / tx.sum()
    # composition control: intra pair distributions re-mixed with the inter pair composition
    fx = np.array([W[1][k].sum() for k in iu])
    fx = fx / fx.sum()
    adj = sum(fx[n] * (W[0][k] / W[0][k].sum()) for n, k in enumerate(iu))
    per = []
    for n, k in enumerate(iu):
        a, b = W[0][k] / W[0][k].sum(), W[1][k] / W[1][k].sum()
        nmin = min(raw[0][k].sum(), raw[1][k].sum())
        pooled = (raw[0][k] + raw[1][k]) / (raw[0][k] + raw[1][k]).sum()
        null = np.mean(
            [
                jsd(
                    rng.multinomial(int(nmin), pooled) / nmin,
                    rng.multinomial(int(nmin), pooled) / nmin,
                )
                for _ in range(30)
            ]
        )
        per.append(
            {
                "jsd": round(jsd(a, b), 5),
                "null": round(float(null), 6),
                "mean_diff": round(float((b * centres).sum() - (a * centres).sum()), 3),
            }
        )
    out["cross"][mode] = {
        "marg_intra": [round(float(x), 5) for x in mi],
        "marg_inter": [round(float(x), 5) for x in mx],
        "adj_intra": [round(float(x), 5) for x in adj],
        "jsd_marg": round(jsd(mi, mx), 5),
        "jsd_composition_adjusted": round(jsd(adj, mx), 5),
        "per_pair": per,
        "mean_intra": round(float((mi * centres).sum()), 3),
        "mean_inter": round(float((mx * centres).sum()), 3),
    }
json.dump(out, open(sys.argv[2], "w"), separators=(",", ":"))
print("written", sys.argv[2])
