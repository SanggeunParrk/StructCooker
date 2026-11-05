import numpy as np
import pytest
from fake_cifmol import backbone_residue, chain, make_mol, protein_chain, residue

from structcooker.instructions.potential import pairs
from structcooker.instructions.potential.pairs import (
    BIN,
    NBIN,
    NV,
    TYPE_INDEX,
    T,
    count_assembly,
    count_assembly_all,
    count_frame_contacts,
    representative_atoms,
    residue_types,
)


def _two_chains(same_cluster=False, n=8):
    a = protein_chain("A", n, cluster="c1")
    b = protein_chain("B", n, z=4.0, sign=-1.0, cluster="c1" if same_cluster else "c2")
    return make_mol([a, b])


def _brute_inter_pairs(mol, cutoff=5.0):
    """Residue pairs of different chains with any two atoms closer than ``cutoff``."""
    xyz = mol.atoms.xyz.value
    res = mol.index_table.atom_res
    chain_of = mol.index_table.res_chain
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    out = set()
    for i, j in zip(*np.nonzero(d < cutoff), strict=True):
        if chain_of[res[i]] < chain_of[res[j]]:
            out.add((int(res[i]), int(res[j])))
    return out


def test_residue_types_for_protein_dna_rna_and_hybrid():
    def nuc(comp, one):
        return residue(comp, one, [("C1'", "C", (0, 0, 0))])

    mol = make_mol(
        [
            chain(
                "P",
                [
                    residue("ALA", "A", []),
                    residue("MSE", "X", []),
                    residue("XXX", "?", []),
                ],
            ),
            chain(
                "D", [nuc("DA", "A"), nuc("DX", "N")], etype="polydeoxyribonucleotide"
            ),
            chain("R", [nuc("U", "U"), nuc("A", "A")], etype="polyribonucleotide"),
            chain(
                "H",
                [nuc("DT", "T"), nuc("U", "U")],
                etype="polydeoxyribonucleotide/polyribonucleotide hybrid",
            ),
            chain("L", [residue("XYZ", "X", [])], etype="non-polymer"),
        ]
    )
    rtype, poly = residue_types(mol)
    names = [pairs.TYPES[t] if t >= 0 else None for t in rtype]
    assert names == ["A", "X", "X", "DA", "DN", "RU", "RA", "DT", "RU", None]
    assert poly.tolist() == [True] * 9 + [False]


def test_representative_atom_is_cb_ca_for_glycine_c4_c2_for_nucleotides():
    prot = chain(
        "P",
        [backbone_residue(0, "A", "ALA"), backbone_residue(1, "G", "GLY", cb=False)],
    )
    na = chain(
        "D",
        [
            residue("DA", "A", [("C4", "C", (0, 0, 0)), ("C2", "C", (1, 0, 0))]),
            residue("DC", "C", [("C4", "C", (2, 0, 0)), ("C2", "C", (3, 0, 0))]),
        ],
        etype="polydeoxyribonucleotide",
    )
    mol = make_mol([prot, na])
    rtype, _ = residue_types(mol)
    rep = representative_atoms(mol, rtype)
    names = mol.atoms.id.value[rep]
    assert names.tolist() == ["CB", "CA", "C4", "C2"]


def test_residue_without_representative_atom_gets_minus_one():
    a = protein_chain("A", 2)
    a["residues"][0]["atoms"] = [x for x in a["residues"][0]["atoms"] if x[0] != "CB"]
    mol = make_mol([a])
    rtype, _ = residue_types(mol)
    assert representative_atoms(mol, rtype)[0] == -1


def test_count_assembly_counts_each_touching_inter_chain_pair_once():
    mol = _two_chains()
    hist = count_assembly(mol)
    assert hist.shape == (2, T, T, NBIN + 1)
    expected = _brute_inter_pairs(mol)
    assert hist[0].sum() == 0
    assert hist[1].sum() == len(expected)
    ala = TYPE_INDEX["A"]
    # residue i of A sits right above residue i of B: CB-CB is 1.0 A -> bin 2
    assert hist[1, ala, ala, 2] == 8
    assert np.flatnonzero(hist[1].sum(axis=(0, 1)))[0] == 2


def test_count_assembly_dedup_collapses_homo_oligomer_copies():
    mol = _two_chains(same_cluster=True)
    n_all = count_assembly(mol, dedup=False)[1].sum()
    n_dedup = count_assembly(mol, dedup=True)[1].sum()
    brute = _brute_inter_pairs(mol)
    positions = {tuple(sorted((a % 8, b % 8))) for a, b in brute}
    assert n_all == len(brute)
    assert n_dedup == len(positions) < n_all


def test_count_assembly_weight_and_accumulation_into_given_hist():
    mol = _two_chains()
    h1 = count_assembly(mol, weight=0.5)
    h2 = count_assembly(mol, weight=0.5, hist=h1)
    assert h2 is h1
    assert h1.sum() == pytest.approx(count_assembly(mol).sum())


def test_count_assembly_min_sep_filters_intra_chain_pairs_by_sequence_distance():
    mol = make_mol([protein_chain("A", 6)])
    xyz, res = mol.atoms.xyz.value, mol.index_table.atom_res
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    touching = {
        (int(res[i]), int(res[j]))
        for i, j in zip(*np.nonzero(d < 5.0), strict=True)
        if res[i] < res[j]
    }
    assert count_assembly(mol, min_sep=3).sum() == sum(b - a >= 3 for a, b in touching)
    h = count_assembly(mol, min_sep=1)
    assert h[1].sum() == 0
    assert h[0].sum() == len(touching) > 0
    assert count_assembly(mol, min_sep=2).sum() == sum(b - a >= 2 for a, b in touching)


def test_count_assembly_without_polymer_returns_zero_hist():
    lig = chain(
        "L", [residue("XYZ", "X", [("C1", "C", (0, 0, 0))])], etype="non-polymer"
    )
    assert count_assembly(make_mol([lig])).sum() == 0
    assert count_assembly_all(make_mol([lig])).sum() == 0


def test_count_assembly_all_uses_rep_atom_distance_within_rmax():
    mol = _two_chains()
    hist = count_assembly_all(mol, rmax=15.0)
    assert hist.shape == (2, T, T, int(15.0 / BIN) + 1)
    ala = TYPE_INDEX["A"]
    # intra: CB(i)-CB(i+3) = 11.4 A -> bin 22, 5 per chain; i+4 is 15.2 A, out of range
    assert hist[0, ala, ala, 22] == 10
    assert hist[0].sum() == 10
    # inter: all (i, j) with sqrt((3.8 d)^2 + 1) < 15
    xyz = mol.atoms.xyz.value
    cb = xyz[mol.atoms.id.value == "CB"]
    d = np.linalg.norm(cb[:8, None] - cb[None, 8:], axis=-1)
    assert hist[1].sum() == (d < 15).sum()
    assert hist[1, ala, ala, 2] == 8
    assert hist[:, :, :, -1].sum() == 0


def test_count_assembly_all_dedup_and_single_residue_cases():
    same = _two_chains(same_cluster=True)
    cb = same.atoms.xyz.value[same.atoms.id.value == "CB"]
    chain_of = np.repeat([0, 1], 8)
    pos = np.tile(np.arange(8), 2)
    d = np.linalg.norm(cb[:, None] - cb[None], axis=-1)
    allowed = (d < 15.0) & (np.arange(16)[:, None] < np.arange(16)[None])
    allowed &= (chain_of[:, None] != chain_of[None]) | (
        abs(pos[:, None] - pos[None]) >= 3
    )
    ii, jj = np.nonzero(allowed)
    unique_positions = {
        tuple(sorted((pos[i], pos[j]))) for i, j in zip(ii, jj, strict=True)
    }
    assert count_assembly_all(same, dedup=False).sum() == len(ii)
    assert count_assembly_all(same, dedup=True).sum() == len(unique_positions)
    assert count_assembly_all(make_mol([protein_chain("A", 1)])).sum() == 0
    far = make_mol([protein_chain("A", 2), protein_chain("B", 2, z=100.0)])
    assert count_assembly_all(far, rmax=5.0)[1].sum() == 0


def _gly_ala():
    gly = chain("A", [backbone_residue(0, "G", "GLY", cb=False)])
    ala = chain("B", [backbone_residue(0, "A", "ALA", z=4.0, sign=-1.0)])
    return make_mol([gly, ala])


def test_frame_contacts_place_the_partner_in_the_reference_frame():
    grid, added, outside = count_frame_contacts(_gly_ala())
    assert grid.shape == (2, T, T, NV, NV, NV)
    assert (added, outside) == (2, 0)
    g, a = TYPE_INDEX["G"], TYPE_INDEX["A"]
    # Gly frame at (0,0,0): partner CB(B) = (0,-.5,2.5) -> local (-0.18, -0.47, 2.5) -> cell (13, 13, 16)
    assert grid[1, g, a, 13, 13, 16] == 1
    # Ala frame at (0,0,4): partner CA(Gly) -> local (0, 0, -4) -> cell (14, 14, 10)
    assert grid[1, a, g, 14, 14, 10] == 1


def test_frame_contacts_weight_goes_to_wgrid_only():
    wgrid = np.zeros((2, T, T, NV, NV, NV), dtype=np.float32)
    grid, _, _ = count_frame_contacts(_gly_ala(), wgrid=wgrid, weight=0.25)
    g, a = TYPE_INDEX["G"], TYPE_INDEX["A"]
    assert grid[1, g, a, 13, 13, 16] == 1
    assert wgrid[1, g, a, 13, 13, 16] == pytest.approx(0.25)


def test_frame_contacts_sidechain_mode_ignores_glycine_and_backbone_only_touch():
    grid, added, outside = count_frame_contacts(_gly_ala(), sidechain=True)
    assert grid.shape == (2, 20, 20, NV, NV, NV)
    assert (added, outside) == (0, 0)


def test_frame_contacts_sidechain_mode_counts_cb_cb_touch_in_both_directions():
    mol = make_mol(
        [
            chain("A", [backbone_residue(0, "A", "ALA")]),
            chain("B", [backbone_residue(0, "A", "ALA", z=4.0, sign=-1.0)]),
        ]
    )
    grid, added, outside = count_frame_contacts(mol, sidechain=True)
    assert (added, outside) == (2, 0)
    assert grid[1, 0, 0, 13, 13, 16] == 1
    assert grid[1, 0, 0, 13, 13, 11] == 1


def test_frame_contacts_count_partners_outside_the_grid_box():
    a = chain("A", [backbone_residue(0, "A", "ALA")])
    far = backbone_residue(0, "A", "ALA", z=4.0, sign=-1.0)
    far["atoms"] = [
        x if x[0] != "CB" else ("CB", "C", (0.0, -0.5, 30.0)) for x in far["atoms"]
    ]
    far["atoms"].append(("CG", "C", (0.0, -0.5, 2.5)))  # touches CB of A
    mol = make_mol([a, chain("B", [far])])
    _, added, outside = count_frame_contacts(mol, sidechain=True)
    assert (added, outside) == (1, 1)


def test_frame_contacts_without_contacts_leave_the_grid_empty():
    mol = make_mol([protein_chain("A", 1), protein_chain("B", 1, z=50.0)])
    grid, added, outside = count_frame_contacts(mol)
    assert (added, outside) == (0, 0)
    assert grid.dtype == np.float32
