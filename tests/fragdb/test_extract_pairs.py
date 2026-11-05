import numpy as np
from fake_cifmol import backbone_residue, chain, make_mol, protein_chain, residue

from structcooker.instructions.fragdb.extract import (
    ATOMS,
    cluster_number,
    extract_pairs,
    residue_coords,
    valid_windows,
    virtual_cb,
)
from structcooker.instructions.mols.backbone import residue_types


def _two(same_cluster=False, n=8):
    a = protein_chain("A", n, cluster="cluster_1")
    b = protein_chain(
        "B", n, z=4.0, sign=-1.0, cluster="cluster_1" if same_cluster else "cluster_2"
    )
    return make_mol([a, b])


def test_cluster_number_takes_the_trailing_integer():
    assert cluster_number("cPDB_000123") == 123
    assert cluster_number("abc") == 0
    assert cluster_number("") == 0
    assert cluster_number("12_x") == 0


def test_residue_coords_have_n_ca_c_o_cb_and_a_virtual_cb_for_glycine():
    ch = chain(
        "A",
        [backbone_residue(0, "A", "ALA"), backbone_residue(1, "G", "GLY", cb=False)],
    )
    mol = make_mol([ch])
    rtype, _ = residue_types(mol)
    x = residue_coords(mol, rtype)
    assert x.shape == (2, len(ATOMS), 3)
    assert np.allclose(x[0, 4], [0.0, -0.5, 1.5])  # the real CB is kept
    assert np.allclose(x[1, 4], virtual_cb(x[1:, 0], x[1:, 1], x[1:, 2])[0])
    assert np.isfinite(x).all()


def test_residue_coords_are_nan_for_non_polymer_residues():
    lig = chain(
        "L", [residue("XYZ", "X", [("N", "N", (0, 0, 0))])], etype="non-polymer"
    )
    mol = make_mol([protein_chain("A", 1), lig])
    rtype, _ = residue_types(mol)
    x = residue_coords(mol, rtype)
    assert np.isfinite(x[0, :4]).all()
    assert np.isnan(x[1, :4]).all()


def test_valid_windows_need_a_complete_contiguous_k_mer():
    mol = make_mol([protein_chain("A", 7)])
    rtype, _ = residue_types(mol)
    x = residue_coords(mol, rtype)
    assert valid_windows(mol, rtype, x, 3).tolist() == [False] + [True] * 5 + [False]
    assert (
        valid_windows(mol, rtype, x, 7).tolist() == [False] * 3 + [True] + [False] * 3
    )
    assert not valid_windows(mol, rtype, x, 9).any()


def test_valid_windows_break_at_a_missing_atom_and_at_a_long_peptide_bond():
    a = protein_chain("A", 7)
    a["residues"][2]["atoms"] = [x for x in a["residues"][2]["atoms"] if x[0] != "O"]
    mol = make_mol([a])
    rtype, _ = residue_types(mol)
    ok = valid_windows(mol, rtype, residue_coords(mol, rtype), 3)
    assert ok.tolist() == [False, False, False, False, True, True, False]

    b = protein_chain("A", 7, x0=0.0)
    for r in b["residues"][4:]:  # shift the C-terminal part: C(3)-N(4) > 2 A
        r["atoms"] = [(n, e, (p[0] + 1.5, p[1], p[2])) for n, e, p in r["atoms"]]
    mol = make_mol([b])
    rtype, _ = residue_types(mol)
    ok = valid_windows(mol, rtype, residue_coords(mol, rtype), 3)
    assert ok.tolist() == [False, True, True, False, False, True, False]


def test_valid_windows_break_at_a_numbering_gap():
    a = protein_chain("A", 7)
    mol = make_mol([a])
    mol.residues.cif_idx.value[4:] += 1
    rtype, _ = residue_types(mol)
    ok = valid_windows(mol, rtype, residue_coords(mol, rtype), 3)
    assert ok.tolist() == [False, True, True, False, False, True, False]


def test_extract_pairs_returns_none_without_a_valid_window_or_contact():
    assert extract_pairs(make_mol([protein_chain("A", 2)]), k=3) is None
    apart = make_mol([protein_chain("A", 6), protein_chain("B", 6, z=60.0)])
    assert extract_pairs(apart, k=3) is None


def test_neighbours_in_one_chain_overlap_and_are_not_pairs():
    assert extract_pairs(make_mol([protein_chain("A", 10)]), k=3) is None


def test_extract_pairs_between_two_chains():
    out = extract_pairs(_two(), k=3)
    n = len(out["coords"])
    # centres 1..6 of each chain; CB(i)-CB'(j) < 5 A only for |i - j| <= 1: 6 + 2 * 5 pairs
    assert n == 16
    assert out["coords"].shape == (n, 2, 3, 5, 3)
    assert out["coords"].dtype == np.int16
    assert out["seq"].shape == (n, 2, 3)
    assert (out["seq"] == 0).all()
    assert out["inter"].tolist() == [1] * n
    assert set(out["chain"].ravel().tolist()) == {"A", "B"}
    assert (out["chain"][:, 0] == "A").all()
    assert (out["cluster"] == [[1, 2]]).all()
    assert out["center"].min() == 1
    assert out["center"].max() == 6
    assert (np.abs(out["center"][:, 0] - out["center"][:, 1]) <= 1).all()
    assert (out["n_atom_pairs"] == 1).all()


def test_coordinates_are_in_the_frame_of_the_first_centre_in_hundredths_of_an_angstrom():
    out = extract_pairs(_two(), k=3)
    ca = ATOMS.index("CA")
    c = ATOMS.index("C")
    # fragment A, centre residue: CA at the origin, C on the +x axis at |C - CA|
    assert (out["coords"][:, 0, 1, ca] == 0).all()
    assert (out["coords"][:, 0, 1, c] == [round(100 * np.hypot(1.3, 0.5)), 0, 0]).all()
    # same-position pair: B's centre CA is 4 A above (frame z axis)
    same = np.flatnonzero(out["center"][:, 0] == out["center"][:, 1])
    assert len(same) == 6
    assert (out["coords"][same, 1, 1, ca] == [0, 0, 400]).all()


def test_homo_oligomer_copies_are_stored_once():
    distinct = extract_pairs(_two(), k=3)
    homo = extract_pairs(_two(same_cluster=True), k=3)
    # (i, j) and (j, i) are the same pair by (cluster, position): 6 + 5 unordered position pairs
    assert len(distinct["coords"]) == 16
    assert len(homo["coords"]) == 11


def test_k_controls_fragment_length_and_centre_range():
    out = extract_pairs(_two(), k=5)
    assert out["coords"].shape[1:] == (2, 5, 5, 3)
    assert out["center"].min() == 2
    assert out["center"].max() == 5
