"""Residue types and representative atoms of a CIFMol.

The representative atom of a residue is AlphaFold3's pseudo-beta (alphafold3/model/protein_data_processing.py): protein CB (CA for glycine); nucleic acids C4 for purines (A, G, DA, DG) and C2 for pyrimidines (C, U, DC, DT).
Residue types, 31 = 28 + 3 unknowns: 20 amino acids, protein X, DNA A C G T + N, RNA A C G U + N; modified residues map to their parent through the CCD canonical
one-letter code.
"""

from __future__ import annotations

from typing import Any

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
TYPES = [*AA, "X", "DA", "DC", "DG", "DT", "DN", "RA", "RC", "RG", "RU", "RN"]
TYPE_INDEX = {t: i for i, t in enumerate(TYPES)}
T = len(TYPES)
_DNA_COMP = {"DA", "DC", "DG", "DT", "DU", "DI", "DN"}


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
