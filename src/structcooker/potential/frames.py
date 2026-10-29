"""Backbone frames of residues and the positions of contacting residues in them.

A residue's frame follows the rigid-group convention of AlphaFold2 / AlphaFold3 (``rigid_from_3_points``):
the first atom gives the x axis, the second is the origin, the third lies in the xy plane.
  protein  (C, CA, N)       origin CA, x towards C, N in the xy plane
  nucleic  (C1', C3', C4')  origin C3', x towards C1', C4' in the xy plane
Positions are expressed in a frame as ``R^T (x - origin)``, so every residue has the same orientation convention
and the contacts of many residues of one type can be pooled in 3D.
"""

from __future__ import annotations

import numpy as np

from structcooker.potential.pairs import AA, TYPES


def frame_atom_names(type_idx: int) -> tuple[str, str, str]:
    """(x-axis atom, origin atom, xy-plane atom) of a residue type."""
    return ("C", "CA", "N") if type_idx < len(AA) + 1 else ("C1'", "C3'", "C4'")


def frames(cm, rtype: np.ndarray):
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
        hit = np.flatnonzero((name == want[atom_res]) & np.isfinite(xyz).all(axis=1) & (rtype[atom_res] >= 0))
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
    ok = np.isfinite(rot).all(axis=(1, 2)) & np.isfinite(o).all(axis=1)
    return o, rot, ok


def to_local(origin: np.ndarray, rot: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Coordinates of ``x`` (n, 3) in the frames (origin (n, 3), rot (n, 3, 3))."""
    return np.einsum("nji,nj->ni", rot, x - origin)
