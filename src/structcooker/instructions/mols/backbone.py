"""Residue types and backbone frames of a CIFMol, shared by the fragment DB and the design set.

The representative atom of a residue is AlphaFold3's pseudo-beta (alphafold3/model/protein_data_processing.py): protein CB (CA for glycine); nucleic acids C4 for purines (A, G, DA, DG) and C2 for pyrimidines (C, U, DC, DT).
Residue types, 31 = 28 + 3 unknowns: 20 amino acids, protein X, DNA A C G T + N, RNA A C G U + N; modified residues map to their parent through the CCD canonical
one-letter code. A residue's frame follows the rigid-group convention of AlphaFold2 / AlphaFold3 (``rigid_from_3_points``): the first atom gives the x axis, the
second is the origin, the third lies in the xy plane.

    protein  (C, CA, N)       origin CA, x towards C, N in the xy plane
    nucleic  (C1', C3', C4')  origin C3', x towards C1', C4' in the xy plane

Positions are expressed in a frame as ``R^T (x - origin)``, so every residue has the same orientation convention.
"""

from __future__ import annotations

from typing import Any

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
TYPES = [*AA, "X", "DA", "DC", "DG", "DT", "DN", "RA", "RC", "RG", "RU", "RN"]
TYPE_INDEX = {t: i for i, t in enumerate(TYPES)}
T = len(TYPES)
_DNA_COMP = {"DA", "DC", "DG", "DT", "DU", "DI", "DN"}

BACKBONE = np.array(["N", "CA", "C", "O", "OXT"])


def residue_types(cm: Any) -> tuple[np.ndarray, np.ndarray]:
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


def representative_atoms(cm: Any, rtype: np.ndarray) -> np.ndarray:
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


def frame_atom_names(type_idx: int) -> tuple[str, str, str]:
    """(x-axis atom, origin atom, xy-plane atom) of a residue type."""
    return ("C", "CA", "N") if type_idx < len(AA) + 1 else ("C1'", "C3'", "C4'")


def frames(cm: Any, rtype: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per residue: origin (R, 3), rotation (R, 3, 3) with columns e1 e2 e3, and a validity mask."""
    at, ix = cm.atoms, cm.index_table
    name, xyz = at.id.value, at.xyz.value
    atom_res = ix.atoms_to_residues(np.arange(len(name)))
    n = len(rtype)
    pts = np.full((3, n, 3), np.nan)
    for slot in range(3):
        want = np.full(n, "", dtype="U3")
        for t in np.unique(rtype[rtype >= 0]):
            want[rtype == t] = frame_atom_names(int(t))[slot]
        hit = np.flatnonzero(
            (name == want[atom_res])
            & np.isfinite(xyz).all(axis=1)
            & (rtype[atom_res] >= 0)
        )
        r, first = np.unique(atom_res[hit], return_index=True)
        pts[slot, r] = xyz[hit[first]]
    a, o, c = pts
    e1 = a - o
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    v = c - o
    e2 = v - (v * e1).sum(axis=1, keepdims=True) * e1
    e2 /= np.linalg.norm(e2, axis=1, keepdims=True)
    e3 = np.cross(e1, e2)
    rot = np.stack([e1, e2, e3], axis=2)
    ok = np.asarray(np.isfinite(rot).all(axis=(1, 2)) & np.isfinite(o).all(axis=1))
    return o, rot, ok


def to_local(origin: np.ndarray, rot: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Coordinates of ``x`` (n, 3) in the frames (origin (n, 3), rot (n, 3, 3))."""
    return np.einsum("nji,nj->ni", rot, x - origin)
