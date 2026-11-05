"""Backbone PDB export of a BioMol item (input of foldseek / mkdssp)."""

from __future__ import annotations

import string
from pathlib import Path

import numpy as np

from structcooker.instructions.design.ss import backbone


CHAIN_LETTERS = string.ascii_uppercase + string.ascii_lowercase + string.digits


def n_protein_chains(cm):
    """Number of chains that have at least one standard amino acid residue with a complete backbone."""
    x, chain, _, rtype, _ = backbone(cm)
    use = (rtype >= 0) & (rtype < 20) & np.isfinite(x).all(axis=(1, 2))
    return len(np.unique(chain[use]))


def backbone_pdb_text(cm):
    """N CA C O of every standard amino acid residue with complete backbone, as PDB text.

    The chains that have such residues get the letters A-Z a-z 0-9 in order (62 at most; more raises ValueError).
    """
    x, chain, cif, rtype, comp = backbone(cm)
    use = (rtype >= 0) & (rtype < 20) & np.isfinite(x).all(axis=(1, 2))
    present, letter = np.unique(chain[use], return_inverse=True)
    if len(present) > len(CHAIN_LETTERS):
        raise ValueError(f"{len(present)} protein chains do not fit the {len(CHAIN_LETTERS)} PDB chain letters")
    chain_letter = np.zeros(len(chain), dtype=np.int64)
    chain_letter[use] = letter
    lines, serial = [], 1
    for r in np.flatnonzero(use):
        for k, a in enumerate(("N", "CA", "C", "O")):
            X, Y, Z = x[r, k]
            lines.append(f"ATOM  {serial % 100000:5d} {a:^4s} {comp[r]:>3s} {CHAIN_LETTERS[chain_letter[r]]}{(cif[r] + 1) % 10000:4d}    "
                         f"{X:8.3f}{Y:8.3f}{Z:8.3f}  1.00  0.00          {a[0]:>2s}")
            serial += 1
    return "\n".join(lines) + "\nEND\n"


def write_backbone_pdb(cm, path):
    """Write ``backbone_pdb_text``; returns the indices of the residues written."""
    Path(path).write_text(backbone_pdb_text(cm))
    x, _, _, rtype, _ = backbone(cm)
    return np.flatnonzero((rtype >= 0) & (rtype < 20) & np.isfinite(x).all(axis=(1, 2)))
