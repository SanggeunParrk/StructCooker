"""Residue-pair distance statistics for a statistical distance potential.

Two ways to choose the pairs (``count_assembly`` = definition A, ``count_assembly_all`` = definition B):

A. a pair is two residues whose heavy atoms touch (within ``contact`` = 5 A): the distribution of the
   representative-atom distance *given contact*, which shows geometry but not how often types touch.
B. a pair is any two residues whose representative atoms are within ``rmax`` (15 A), touching or not (the
   DOPE / DFIRE convention): the distance histogram then also carries how often a type pair comes close.

For every pair the distance between the two *representative atoms* is binned, by the unordered residue-type
pair and by scope: ``intra`` (same chain, sequence separation >= ``min_sep``) or ``inter`` (different chains).

Representative atom = AlphaFold3's pseudo-beta (alphafold3/model/protein_data_processing.py): protein CB
(CA for glycine); nucleic acids C4 for purines (A, G, DA, DG) and C2 for pyrimidines (C, U, DC, DT).
Residue types, 31 = 28 + 3 unknowns: 20 amino acids, protein X, DNA A C G T + N, RNA A C G U + N.
Modified residues map to their parent through the CCD canonical one-letter code. A residue without a
representative atom (missing / no coordinates) is skipped.

Within one assembly a pair that is identical by (30% cluster, residue position) on both sides is counted once, so
the copies of a homo-oligomer (chains of one cluster, or the symmetry-equivalent interfaces) do not multiply
a pair.
"""

from __future__ import annotations

import numpy as np

from structcooker.instructions.interface.detect import inter_node_pairs

AA = "ACDEFGHIKLMNPQRSTVWY"
TYPES = [*AA, "X", "DA", "DC", "DG", "DT", "DN", "RA", "RC", "RG", "RU", "RN"]
TYPE_INDEX = {t: i for i, t in enumerate(TYPES)}
T = len(TYPES)
BIN, RMAX = 0.5, 20.0
NBIN = int(RMAX / BIN)          # the last bin (index NBIN) is the overflow >= RMAX
SCOPES = ("intra", "inter")
_DNA_COMP = {"DA", "DC", "DG", "DT", "DU", "DI", "DN"}


def residue_types(cm) -> tuple[np.ndarray, np.ndarray]:
    """``(type index (R,), is_polymer (R,))``; -1 for residues that are not polymer residues."""
    rs, ch, ix = cm.residues, cm.chains, cm.index_table
    one = rs.one_letter_code_can.value
    comp = rs.chem_comp_id.value
    chain = ix.residues_to_chains(np.arange(len(one)))
    etype = ch.entity_type.value[chain]
    out = np.full(len(one), -1, dtype=np.int16)
    for r in range(len(one)):
        et, letter = str(etype[r]), str(one[r])
        if et.startswith("polypeptide"):
            out[r] = TYPE_INDEX[letter] if letter in AA else TYPE_INDEX["X"]
        elif et.startswith(("polydeoxyribo", "polyribo")):
            dna = et.startswith("polydeoxyribo")
            if "hybrid" in et:  # DNA / RNA is decided per residue
                dna = str(comp[r]) in _DNA_COMP
            letters = "ACGT" if dna else "ACGU"
            code = letter if letter in letters else "N"
            out[r] = TYPE_INDEX[("D" if dna else "R") + code]
    return out, out >= 0


def representative_atoms(cm, rtype: np.ndarray) -> np.ndarray:
    """Atom index of each residue's pseudo-beta atom, -1 if absent or without coordinates."""
    at, ix = cm.atoms, cm.index_table
    name = at.id.value
    xyz = at.xyz.value
    atom_res = ix.atoms_to_residues(np.arange(len(name)))
    target = np.full(len(rtype), "", dtype="U3")
    for r, t in enumerate(rtype.tolist()):
        if t < 0:
            continue
        label = TYPES[t]
        if t < len(AA) + 1:
            target[r] = "CA" if label == "G" else "CB"
        else:
            target[r] = "C4" if label[1] in ("A", "G") else "C2"
    hit = np.flatnonzero((name == target[atom_res]) & np.isfinite(xyz).all(axis=1))
    rep = np.full(len(rtype), -1, dtype=np.int64)
    r_hit, first = np.unique(atom_res[hit], return_index=True)
    rep[r_hit] = hit[first]
    return rep


def count_assembly(cm, contact: float = 5.0, min_sep: int = 3, weight: float = 1.0, hist=None, dedup: bool = True):
    """Add this assembly's pairs to ``hist`` (``(2, T, T, NBIN + 1)`` float64, created if None)."""
    if hist is None:
        hist = np.zeros((2, T, T, NBIN + 1))
    rs, ch, at, ix = cm.residues, cm.chains, cm.atoms, cm.index_table
    rtype, poly = residue_types(cm)
    if not poly.any():
        return hist
    rep = representative_atoms(cm, rtype)
    xyz = at.xyz.value
    atom_res = ix.atoms_to_residues(np.arange(len(xyz)))
    use = poly[atom_res] & np.isfinite(xyz).all(axis=1)
    idx = np.flatnonzero(use)
    i, j, _ = inter_node_pairs(xyz[idx], atom_res[idx], contact)
    if len(i) == 0:
        return hist
    ra, rb = atom_res[idx[i]], atom_res[idx[j]]
    lo, hi = np.minimum(ra, rb), np.maximum(ra, rb)
    nres = len(rtype)
    key = np.unique(lo.astype(np.int64) * nres + hi)
    lo, hi = key // nres, key % nres
    ok = (rep[lo] >= 0) & (rep[hi] >= 0)
    lo, hi = lo[ok], hi[ok]
    chain = ix.residues_to_chains(np.arange(nres))
    cif = rs.cif_idx.value.astype(np.int64)
    same = chain[lo] == chain[hi]
    keep = ~same | (np.abs(cif[lo] - cif[hi]) >= min_sep)
    lo, hi, same = lo[keep], hi[keep], same[keep]
    if dedup and len(lo):
        cl = ch.cluster_id.value
        _, cint = np.unique(cl, return_inverse=True)
        cint = cint.astype(np.int64)
        ca, cb = cint[chain[lo]], cint[chain[hi]]
        a = ca * 100000 + cif[lo]
        b = cb * 100000 + cif[hi]
        first = np.minimum(a, b)
        second = np.maximum(a, b)
        _, u = np.unique(np.stack([first, second], axis=1), axis=0, return_index=True)
        lo, hi, same = lo[u], hi[u], same[u]
    d = np.linalg.norm(xyz[rep[lo]] - xyz[rep[hi]], axis=1)
    b = np.minimum((d / BIN).astype(np.int64), NBIN)
    ta, tb = rtype[lo].astype(np.int64), rtype[hi].astype(np.int64)
    np.add.at(hist, (np.where(same, 0, 1), np.minimum(ta, tb), np.maximum(ta, tb), b), weight)
    return hist


def count_assembly_all(cm, rmax: float = 15.0, min_sep: int = 3, weight: float = 1.0, hist=None, dedup: bool = True):
    """Definition B: every residue pair whose representative atoms are within ``rmax``.

    ``hist`` is ``(2, T, T, int(rmax / BIN) + 1)``; the last bin is never filled (distances are < ``rmax``).
    """
    nb = int(rmax / BIN)
    if hist is None:
        hist = np.zeros((2, T, T, nb + 1))
    rs, ch, ix = cm.residues, cm.chains, cm.index_table
    rtype, poly = residue_types(cm)
    if not poly.any():
        return hist
    rep = representative_atoms(cm, rtype)
    have = np.flatnonzero(rep >= 0)
    if len(have) < 2:
        return hist
    xyz = cm.atoms.xyz.value[rep[have]]
    i, j, d = inter_node_pairs(xyz, np.arange(len(have)), rmax)
    if len(i) == 0:
        return hist
    lo, hi = have[i], have[j]                      # have is ascending, so lo < hi
    chain = ix.residues_to_chains(np.arange(len(rtype)))
    cif = rs.cif_idx.value.astype(np.int64)
    same = chain[lo] == chain[hi]
    keep = ~same | (np.abs(cif[lo] - cif[hi]) >= min_sep)
    lo, hi, same, d = lo[keep], hi[keep], same[keep], d[keep]
    if dedup and len(lo):
        _, cint = np.unique(ch.cluster_id.value, return_inverse=True)
        cint = cint.astype(np.int64)
        a = cint[chain[lo]] * 100000 + cif[lo]
        b = cint[chain[hi]] * 100000 + cif[hi]
        key = np.minimum(a, b) * (1 << 28) + np.maximum(a, b)   # one integer per unordered (cluster, position) pair
        _, u = np.unique(key, return_index=True)
        lo, hi, same, d = lo[u], hi[u], same[u], d[u]
    b = np.minimum((d / BIN).astype(np.int64), nb)
    ta, tb = rtype[lo].astype(np.int64), rtype[hi].astype(np.int64)
    np.add.at(hist, (np.where(same, 0, 1), np.minimum(ta, tb), np.maximum(ta, tb), b), weight)
    return hist


# ---------------------------------------------------------------------------------------------------------------
# contacts seen from each residue's own backbone frame

BOX, VOX = 14.0, 1.0                    # the grid spans +-BOX A around the frame origin
NV = int(2 * BOX / VOX)                 # voxels per axis


BACKBONE = np.array(["N", "CA", "C", "O", "OXT"])


def count_frame_contacts(cm, contact: float = 5.0, min_sep: int = 3, grid=None, wgrid=None, weight: float = 1.0,
                         sidechain: bool = False):
    """For every touching residue pair (definition A), add the partner's representative atom to the grid of the
    other residue's frame, in both directions. ``grid`` is ``(2, T, T, NV, NV, NV)`` float32 indexed by
    ``[scope, reference type, partner type, x, y, z]``; ``wgrid`` (same shape) receives the same samples times
    ``weight``. Returns ``(grid, n_added, n_outside)``.

    ``sidechain=True`` narrows the contact: only the 20 standard amino acids count and two residues touch when a
    side-chain heavy atom of one (every heavy atom except N, CA, C, O, OXT; Gly has none) is within ``contact`` of a
    side-chain heavy atom of the other. Then ``grid`` is ``(2, 20, 20, NV, NV, NV)``.
    """
    from structcooker.instructions.potential.frames import frames, to_local

    nt = len(AA) if sidechain else T
    if grid is None:
        grid = np.zeros((2, nt, nt, NV, NV, NV), dtype=np.float32)
    rs, ch, at, ix = cm.residues, cm.chains, cm.atoms, cm.index_table
    rtype, poly = residue_types(cm)
    if sidechain:
        poly = poly & (rtype < len(AA))             # the 20 standard amino acids only
    if not poly.any():
        return grid, 0, 0
    rep = representative_atoms(cm, rtype)
    origin, rot, ok = frames(cm, rtype)
    xyz = at.xyz.value
    atom_res = ix.atoms_to_residues(np.arange(len(xyz)))
    elig = poly[atom_res] & np.isfinite(xyz).all(axis=1)
    if sidechain:
        elig &= ~np.isin(at.id.value, BACKBONE)
    use = np.flatnonzero(elig)
    i, j, _ = inter_node_pairs(xyz[use], atom_res[use], contact)
    if len(i) == 0:
        return grid, 0, 0
    ra, rb = atom_res[use[i]], atom_res[use[j]]
    nres = len(rtype)
    key = np.unique(np.minimum(ra, rb).astype(np.int64) * nres + np.maximum(ra, rb))
    lo, hi = key // nres, key % nres
    chain = ix.residues_to_chains(np.arange(nres))
    cif = rs.cif_idx.value.astype(np.int64)
    same = chain[lo] == chain[hi]
    keep = ~same | (np.abs(cif[lo] - cif[hi]) >= min_sep)
    lo, hi, same = lo[keep], hi[keep], same[keep]
    # both directions: (reference, partner)
    ref = np.concatenate([lo, hi])
    par = np.concatenate([hi, lo])
    sc = np.concatenate([same, same])
    good = ok[ref] & (rep[par] >= 0)
    ref, par, sc = ref[good], par[good], sc[good]
    if len(ref) == 0:
        return grid, 0, 0
    _, cint = np.unique(ch.cluster_id.value, return_inverse=True)
    cint = cint.astype(np.int64)
    a = cint[chain[ref]] * 100000 + cif[ref]
    b = cint[chain[par]] * 100000 + cif[par]
    _, u = np.unique(a * (1 << 28) + b, return_index=True)       # symmetry copies: one sample per (cluster, position) pair
    ref, par, sc = ref[u], par[u], sc[u]
    local = to_local(origin[ref], rot[ref], xyz[rep[par]])
    cell = np.floor((local + BOX) / VOX).astype(np.int64)
    inside = ((cell >= 0) & (cell < NV)).all(axis=1)
    idx = (np.where(sc, 0, 1), rtype[ref].astype(np.int64), rtype[par].astype(np.int64), cell[:, 0], cell[:, 1], cell[:, 2])
    sel = tuple(v[inside] for v in idx)
    np.add.at(grid, sel, 1.0)
    if wgrid is not None:
        np.add.at(wgrid, sel, weight)
    return grid, int(inside.sum()), int((~inside).sum())
