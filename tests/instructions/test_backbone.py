"""Residue types and representative atoms of a CIFMol."""

import numpy as np
from fake_cifmol import make_mol, protein_chain

from structcooker.instructions.mols.backbone import representative_atoms, residue_types


def test_representative_atom_is_cb_and_ca_for_glycine():
    mol = make_mol([protein_chain("A", 2, seq="AG")])
    rtype, _ = residue_types(mol)
    rep = representative_atoms(mol, rtype)
    names = mol.atoms.id.value
    assert names[rep[0]] == "CB"
    assert names[rep[1]] == "CA"


def test_residue_without_its_representative_atom_has_none():
    mol = make_mol([protein_chain("A", 1)])
    rtype, _ = residue_types(mol)
    mol.atoms.id.value = np.where(mol.atoms.id.value == "CB", "XX", mol.atoms.id.value)
    assert representative_atoms(mol, rtype)[0] == -1
