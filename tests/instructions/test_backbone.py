import numpy as np
from fake_cifmol import chain, make_mol, protein_chain, residue

from structcooker.instructions.mols.backbone import (
    TYPE_INDEX,
    frame_atom_names,
    frames,
    residue_types,
    to_local,
)


def test_frame_atoms_depend_on_polymer_kind():
    assert frame_atom_names(TYPE_INDEX["A"]) == ("C", "CA", "N")
    assert frame_atom_names(TYPE_INDEX["X"]) == ("C", "CA", "N")
    assert frame_atom_names(TYPE_INDEX["DA"]) == ("C1'", "C3'", "C4'")
    assert frame_atom_names(TYPE_INDEX["RU"]) == ("C1'", "C3'", "C4'")


def test_protein_frame_is_right_handed_with_origin_at_ca():
    mol = make_mol([protein_chain("A", 3)])
    rtype, _ = residue_types(mol)
    origin, rot, ok = frames(mol, rtype)
    assert ok.all()
    ca = mol.atoms.xyz.value[1::5]
    assert np.allclose(origin, ca)
    eye = np.einsum("nji,njk->nik", rot, rot)
    assert np.allclose(eye, np.eye(3), atol=1e-12)
    assert np.allclose(np.linalg.det(rot), 1.0)
    c = mol.atoms.xyz.value[2::5]
    e1 = (c - ca) / np.linalg.norm(c - ca, axis=1, keepdims=True)
    assert np.allclose(rot[:, :, 0], e1)


def test_to_local_puts_c_on_the_x_axis_and_n_in_the_xy_plane():
    mol = make_mol([protein_chain("A", 2)])
    rtype, _ = residue_types(mol)
    origin, rot, _ = frames(mol, rtype)
    xyz = mol.atoms.xyz.value
    c = to_local(origin, rot, xyz[[2, 7]])
    n = to_local(origin, rot, xyz[[0, 5]])
    assert np.allclose(c[:, 1:], 0, atol=1e-12)
    assert np.allclose(c[:, 0], np.hypot(1.3, 0.5))
    assert np.allclose(n[:, 2], 0, atol=1e-12)
    assert (n[:, 1] > 0).all()


def test_nucleic_frame_uses_c1_c3_c4():
    atoms = [("C1'", "C", (1, 0, 0)), ("C3'", "C", (0, 0, 0)), ("C4'", "C", (0, 2, 0))]
    mol = make_mol(
        [chain("D", [residue("DA", "A", atoms)], etype="polydeoxyribonucleotide")]
    )
    rtype, _ = residue_types(mol)
    origin, rot, ok = frames(mol, rtype)
    assert ok.tolist() == [True]
    assert np.allclose(origin, 0)
    assert np.allclose(rot[0], np.eye(3))


def test_residue_missing_a_frame_atom_is_invalid():
    a = protein_chain("A", 2)
    a["residues"][1]["atoms"] = [x for x in a["residues"][1]["atoms"] if x[0] != "N"]
    lig = chain(
        "L", [residue("XYZ", "X", [("C1", "C", (0, 0, 9))])], etype="non-polymer"
    )
    mol = make_mol([a, lig])
    rtype, _ = residue_types(mol)
    _, rot, ok = frames(mol, rtype)
    assert ok.tolist() == [True, False, False]
    assert np.isnan(rot[1]).any()


def test_to_local_inverts_the_frame():
    mol = make_mol([protein_chain("A", 2)])
    rtype, _ = residue_types(mol)
    origin, rot, _ = frames(mol, rtype)
    x = np.array([[1.0, 2.0, 3.0], [-4.0, 0.5, 2.0]])
    local = to_local(origin, rot, x)
    back = np.einsum("nij,nj->ni", rot, local) + origin
    assert np.allclose(back, x)


def test_representative_atom_is_cb_and_ca_for_glycine():
    from structcooker.instructions.mols.backbone import representative_atoms

    mol = make_mol([protein_chain("A", 2, seq="AG")])
    rtype, _ = residue_types(mol)
    rep = representative_atoms(mol, rtype)
    names = mol.atoms.id.value
    assert names[rep[0]] == "CB"
    assert names[rep[1]] == "CA"


def test_residue_without_its_representative_atom_has_none():
    from structcooker.instructions.mols.backbone import representative_atoms

    mol = make_mol([protein_chain("A", 1)])
    rtype, _ = residue_types(mol)
    mol.atoms.id.value = np.where(mol.atoms.id.value == "CB", "XX", mol.atoms.id.value)
    assert representative_atoms(mol, rtype)[0] == -1
