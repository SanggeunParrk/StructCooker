"""Fragment template pairs: two touching k-mers of a protein, with coordinates.

A *fragment* is k consecutive standard amino acids (k = 3, 5, 7) of one chain, all with N, CA, C, O resolved and joined by
peptide bonds (C-N < 2.0 A). A *pair* is two fragments whose *center* residues touch: a side-chain heavy atom of one is within
5 A of a side-chain heavy atom of the other (Gly has no side chain, so it is never a center). Inside one chain the two
fragments must not overlap (centers at least k residues apart in the sequence).

Each pair is stored once, with A the earlier residue (a query for (b, a) is answered from the stored (a, b) by re-expressing the
coordinates in B's frame). Coordinates are in the backbone frame of A's center residue (AF2 convention, see ``potential.frames``),
as int16 in units of 0.01 A: shape ``(2, k, 5, 3)`` for the atoms N, CA, C, O, CB of fragment A then B. CB is the real Cβ, or for
Gly the virtual Cβ built from N, CA, C. Within one assembly a pair that is identical by (30% cluster, center position) on both
sides is stored once (copies of a homo-oligomer).
"""

from __future__ import annotations

import re

import numpy as np

from structcooker.instructions.interface.detect import inter_node_pairs
from structcooker.instructions.potential.frames import frames, to_local
from structcooker.instructions.potential.pairs import AA, BACKBONE, residue_types

ATOMS = ("N", "CA", "C", "O", "CB")
SCALE = 100.0                     # int16 coordinates in 0.01 A
PEPTIDE_BOND_MAX = 2.0


def virtual_cb(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Cβ from the backbone (the AlphaFold2 construction), for glycine and for residues without one."""
    b, cc = ca - n, c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def residue_coords(cm, rtype: np.ndarray) -> np.ndarray:
    """(R, 5, 3) coordinates of N, CA, C, O, CB per residue (NaN where missing; CB virtual when absent)."""
    at, ix = cm.atoms, cm.index_table
    name, xyz = at.id.value, at.xyz.value
    atom_res = ix.atoms_to_residues(np.arange(len(name)))
    out = np.full((len(rtype), len(ATOMS), 3), np.nan)
    fin = np.isfinite(xyz).all(axis=1)
    for k, a in enumerate(ATOMS):
        hit = np.flatnonzero((name == a) & fin & (rtype[atom_res] >= 0))
        r, first = np.unique(atom_res[hit], return_index=True)
        out[r, k] = xyz[hit[first]]
    cb = virtual_cb(out[:, 0], out[:, 1], out[:, 2])
    missing = ~np.isfinite(out[:, 4]).all(axis=1)
    out[missing, 4] = cb[missing]
    return out


def valid_windows(cm, rtype: np.ndarray, x: np.ndarray, k: int) -> np.ndarray:
    """Boolean (R,): a k-mer centred on this residue is complete, contiguous and made of standard amino acids."""
    rs, ix = cm.residues, cm.index_table
    n, half = len(rtype), k // 2
    chain = ix.residues_to_chains(np.arange(n))
    cif = rs.cif_idx.value.astype(np.int64)
    ok_res = (rtype >= 0) & (rtype < len(AA)) & np.isfinite(x[:, :4]).all(axis=(1, 2))
    link = np.zeros(n, dtype=bool)                                   # link[r]: residue r is joined to r + 1
    if n > 1:
        same = (chain[:-1] == chain[1:]) & (cif[1:] - cif[:-1] == 1) & ok_res[:-1] & ok_res[1:]
        d = np.linalg.norm(x[:-1, 2] - x[1:, 0], axis=1)             # C(i) - N(i+1)
        link[:-1] = same & (d < PEPTIDE_BOND_MAX)
    cs_res = np.concatenate([[0], np.cumsum(ok_res)])
    cs_link = np.concatenate([[0], np.cumsum(link)])
    out = np.zeros(n, dtype=bool)
    r = np.arange(half, n - half)
    out[r] = ((cs_res[r + half + 1] - cs_res[r - half]) == k) & ((cs_link[r + half] - cs_link[r - half]) == k - 1)
    return out


def cluster_number(s: str) -> int:
    m = re.search(r"(\d+)$", str(s))
    return int(m.group(1)) if m else 0


def extract_pairs(cm, k: int = 7, contact: float = 5.0):
    """Pairs of touching k-mers of one assembly -> dict of arrays, or None. See the module docstring."""
    rs, ch, at, ix = cm.residues, cm.chains, cm.atoms, cm.index_table
    rtype, poly = residue_types(cm)
    x = residue_coords(cm, rtype)
    ok = valid_windows(cm, rtype, x, k)
    if ok.sum() < 2:
        return None
    origin, rot, fok = frames(cm, rtype)
    ok &= fok
    xyz = at.xyz.value
    atom_res = ix.atoms_to_residues(np.arange(len(xyz)))
    elig = ok[atom_res] & np.isfinite(xyz).all(axis=1) & ~np.isin(at.id.value, BACKBONE)
    use = np.flatnonzero(elig)
    if len(use) < 2:
        return None
    i, j, _ = inter_node_pairs(xyz[use], atom_res[use], contact)
    if len(i) == 0:
        return None
    ra, rb = atom_res[use[i]], atom_res[use[j]]
    nres = len(rtype)
    key, counts = np.unique(np.minimum(ra, rb).astype(np.int64) * nres + np.maximum(ra, rb), return_counts=True)
    lo, hi = key // nres, key % nres
    chain = ix.residues_to_chains(np.arange(nres))
    cif = rs.cif_idx.value.astype(np.int64)
    same = chain[lo] == chain[hi]
    keep = ~same | (np.abs(cif[lo] - cif[hi]) >= k)
    lo, hi, same, counts = lo[keep], hi[keep], same[keep], counts[keep]
    if len(lo) == 0:
        return None
    _, cint = np.unique(ch.cluster_id.value, return_inverse=True)
    cint = cint.astype(np.int64)
    a = cint[chain[lo]] * 100000 + cif[lo]
    b = cint[chain[hi]] * 100000 + cif[hi]
    _, u = np.unique(np.minimum(a, b) * (1 << 28) + np.maximum(a, b), return_index=True)
    lo, hi, same, counts = lo[u], hi[u], same[u], counts[u]
    half = k // 2
    off = np.arange(-half, half + 1)
    wa, wb = lo[:, None] + off, hi[:, None] + off                         # (n, k) residue indices
    xa, xb = x[wa], x[wb]                                                  # (n, k, 5, 3)
    o, r = origin[lo], rot[lo]
    la = np.einsum("nji,nkaj->nkai", r, xa - o[:, None, None, :])
    lb = np.einsum("nji,nkaj->nkai", r, xb - o[:, None, None, :])
    coords = np.clip(np.round(np.stack([la, lb], axis=1) * SCALE), -32767, 32767).astype(np.int16)   # (n, 2, k, 5, 3)
    cl = ch.cluster_id.value
    return {
        "seq": np.stack([rtype[wa], rtype[wb]], axis=1).astype(np.uint8),   # (n, 2, k) amino-acid indices in ``AA`` order
        "coords": coords,
        "center": np.stack([cif[lo], cif[hi]], axis=1).astype(np.int32),   # sequence position of each center
        "chain": np.stack([ch.chain_id.value[chain[lo]], ch.chain_id.value[chain[hi]]], axis=1).astype("U6"),
        "cluster": np.array([[cluster_number(cl[chain[p]]), cluster_number(cl[chain[q]])] for p, q in zip(lo, hi)], dtype=np.uint32).reshape(-1, 2),
        "inter": (~same).astype(np.uint8),
        "n_atom_pairs": np.minimum(counts, 255).astype(np.uint8),
    }
