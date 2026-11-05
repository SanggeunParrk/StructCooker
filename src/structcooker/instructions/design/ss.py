"""Secondary structure of a BioMol item: a vectorised DSSP (Kabsch & Sander 1983) on the backbone atoms.

Hydrogen-bond energy of CO(a) with NH(b): E = 0.084 * 332 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) kcal/mol, a bond if E < -0.5; every
CO and every NH keeps its two best partners. Turns (3, 4, 5), helices (H alpha, G 3-10, I pi), bridges (parallel / antiparallel) and
ladders give the 8 states; ``to_three`` folds them to helix (H G I), strand (E B) and coil (everything else). Residues are neighbours
when they follow each other in one chain and the C-N peptide bond is shorter than 2 A. Not implemented: bulge merging of ladders (a
short bulge inside a sheet can leave one coil residue between two strands) and the turn / bend classes (T S), which are coil in the
three-state reduction anyway. Checked against mkdssp in scripts/design/check_ss.py.
"""

from __future__ import annotations

import numpy as np

from structcooker.instructions.interface.detect import inter_node_pairs
from structcooker.instructions.potential.pairs import residue_types

Q1Q2F = 0.084 * 332.0
HB_CUT = -0.5
PEPTIDE = 2.0
CA_RANGE = 9.0


def backbone(cm):
    """``(x (R, 4, 3) N CA C O or NaN, chain (R,), cif index (R,), residue type (R,), chem_comp_id (R,))``."""
    at, ix, rs = cm.atoms, cm.index_table, cm.residues
    rtype, _ = residue_types(cm)
    name, xyz = at.id.value, at.xyz.value
    atom_res = ix.atoms_to_residues(np.arange(len(name)))
    x = np.full((len(rtype), 4, 3), np.nan)
    fin = np.isfinite(xyz).all(axis=1) & (rtype[atom_res] >= 0)
    for k, a in enumerate(("N", "CA", "C", "O")):
        hit = np.flatnonzero((name == a) & fin)
        r, first = np.unique(atom_res[hit], return_index=True)
        x[r, k] = xyz[hit[first]]
    chain = ix.residues_to_chains(np.arange(len(rtype)))
    return x, chain, rs.cif_idx.value.astype(np.int64), rtype, rs.chem_comp_id.value


def _best2(who, e):
    """Mask of the two lowest-energy entries for every value of ``who``."""
    order = np.lexsort((e, who))
    w = who[order]
    start = np.concatenate([[0], np.flatnonzero(w[1:] != w[:-1]) + 1])
    rank = np.arange(len(w)) - np.repeat(start, np.diff(np.concatenate([start, [len(w)]])))
    keep = np.zeros(len(w), dtype=bool)
    keep[order[rank < 2]] = True
    return keep


def assign8(x, chain, cif, proline, usable):
    """8-state codes (``H G I E B -``) for every residue; unusable residues are ``-``."""
    n = len(chain)
    link = np.zeros(n, dtype=bool)                                      # link[r]: r and r + 1 are joined by a peptide bond
    if n > 1:
        d = np.linalg.norm(x[:-1, 2] - x[1:, 0], axis=1)
        link[:-1] = (chain[:-1] == chain[1:]) & (cif[1:] - cif[:-1] == 1) & usable[:-1] & usable[1:] & (d < PEPTIDE)
    cs_link = np.concatenate([[0], np.cumsum(link)])
    state = np.full(n, "-", dtype="U1")
    ca = x[:, 1]
    ok = np.flatnonzero(usable & np.isfinite(ca).all(axis=1))
    if len(ok) < 5:
        return state
    has_h = np.concatenate([[False], link[:-1]]) & ~proline & usable
    hpos = np.full((n, 3), np.nan)
    v = x[:-1, 2] - x[:-1, 3]                                           # C(i-1) - O(i-1)
    hpos[1:] = x[1:, 0] + v / np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-6)
    hpos[~has_h] = np.nan
    i, j, _ = inter_node_pairs(ca[ok], np.arange(len(ok)), CA_RANGE)
    lo, hi = ok[i], ok[j]                                               # residue pairs with CA-CA < 9 A
    a, b = np.concatenate([lo, hi]), np.concatenate([hi, lo])           # a: acceptor (CO), b: donor (NH)
    keep = (b != a) & (b != a + 1)
    a, b = a[keep], b[keep]
    dist = lambda p, q: np.linalg.norm(p - q, axis=1)                   # noqa: E731
    with np.errstate(divide="ignore", invalid="ignore"):
        e = Q1Q2F * (1 / dist(x[a, 3], x[b, 0]) + 1 / dist(x[a, 2], hpos[b]) - 1 / dist(x[a, 3], hpos[b]) - 1 / dist(x[a, 2], x[b, 0]))
    sel = np.isfinite(e) & (e < HB_CUT)
    a, b, e = a[sel], b[sel], e[sel]
    for who in ("a", "b"):
        k = _best2(a if who == "a" else b, e)
        a, b, e = a[k], b[k], e[k]
    hb = np.sort(a.astype(np.int64) * n + b)                            # hb(a, b): CO(a) ... NH(b)

    def H(p, q):
        p, q = np.asarray(p, dtype=np.int64), np.asarray(q, dtype=np.int64)
        inside = (p >= 0) & (p < n) & (q >= 0) & (q < n)
        key = np.where(inside, p * n + q, -1)
        pos = np.searchsorted(hb, key)
        return inside & (pos < len(hb)) & (hb[np.minimum(pos, max(len(hb) - 1, 0))] == key) if len(hb) else np.zeros(p.shape, dtype=bool)

    def chained(lo_, k):                                                # residues lo_ .. lo_ + k all linked in a row
        lo_ = np.asarray(lo_, dtype=np.int64)
        inside = (lo_ >= 0) & (lo_ + k <= n - 1)
        out = np.zeros(lo_.shape, dtype=bool)
        out[inside] = (cs_link[lo_[inside] + k] - cs_link[lo_[inside]]) == k
        return out

    idx = np.arange(n)
    for m, code in ((4, "H"), (3, "G"), (5, "I")):
        t = H(idx, idx + m) & chained(idx, m)                           # an m-turn starts at i
        for s in np.flatnonzero(t[1:] & t[:-1]) + 1:                    # two consecutive m-turns: residues s .. s + m - 1
            seg = np.arange(s, min(s + m, n))
            state[seg[state[seg] == "-"]] = code
    p, q = np.minimum(lo, hi), np.maximum(lo, hi)
    far = (q - p > 2) | (chain[p] != chain[q])
    p, q = p[far], q[far]
    interior = lambda r: (r >= 1) & (r <= n - 2) & chained(r - 1, 2)    # noqa: E731
    base = interior(p) & interior(q)
    par = base & ((H(p - 1, q) & H(q, p + 1)) | (H(q - 1, p) & H(p, q + 1)))
    anti = base & ((H(p, q) & H(q, p)) | (H(p - 1, q + 1) & H(q - 1, p + 1)))
    bridges = {(int(u), int(w), 0) for u, w in zip(p[par], q[par])} | {(int(u), int(w), 1) for u, w in zip(p[anti], q[anti])}
    in_bridge, in_ladder = np.zeros(n, dtype=bool), np.zeros(n, dtype=bool)
    for u, w, t in bridges:
        in_bridge[[u, w]] = True
        nu, nw = (u + 1, w + 1) if t == 0 else (u + 1, w - 1)
        if (min(nu, nw), max(nu, nw), t) in bridges:                    # the next bridge of the same ladder
            in_ladder[[u, w, nu, nw]] = True
    free = state == "-"
    state[free & in_ladder] = "E"
    state[(state == "-") & in_bridge] = "B"
    state[~usable] = "-"
    return state


def to_three(s8):
    return np.where(np.isin(s8, ["H", "G", "I"]), "H", np.where(np.isin(s8, ["E", "B"]), "E", "C"))


def max_run(mask):
    best = cur = 0
    for v in mask:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


def secondary_structure(cm):
    """``(codes per residue (R,): 'H' helix, 'E' strand, 'C' coil, '-' not assessed, summary dict)`` of a CIFMol / CIFMolAttached.

    Summary: ``helix``, ``strand``, ``coil`` fractions over the residues DSSP could assess (standard amino acids with N, CA, C, O),
    ``max_coil_run`` (longest run of consecutive coil residues), ``n_assessed``.
    """
    x, chain, cif, rtype, comp = backbone(cm)
    usable = (rtype >= 0) & (rtype < 20) & np.isfinite(x).all(axis=(1, 2))
    s3 = to_three(assign8(x, chain, cif, comp == "PRO", usable))
    s3[~usable] = "-"
    n = int(usable.sum())
    frac = lambda c: float((s3 == c).sum() / n) if n else 0.0            # noqa: E731
    return s3, {"helix": frac("H"), "strand": frac("E"), "coil": frac("C"), "max_coil_run": max_run(s3 == "C"), "n_assessed": n}
