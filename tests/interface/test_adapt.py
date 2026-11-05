import numpy as np
from fake_cifmol import chain, make_mol, protein_chain, residue

from structcooker.instructions.interface.adapt import build_assembly, pick_alt_keys

EXCLUDED = frozenset({"HEM_X"})


def _lig(cid="L", comp="XYZ", n=7, z=2.0):
    atoms = [(f"C{i}", "C", (3.8 * i + 1, -2.0, z)) for i in range(n)]
    return chain(cid, [residue(comp, "X", atoms)], etype="non-polymer")


def _single(cid, comp, elem, p):
    return chain(cid, [residue(comp, "X", [(comp, elem, p)])], etype="non-polymer")


def test_pick_alt_keys_prefers_the_no_altloc_entry_per_assembly_model():
    keys = ["1_1_A", "1_1_.", "1_1_B", "1_2_A", "1_2_B", "2_1_B"]
    assert sorted(pick_alt_keys(keys)) == ["1_1_.", "1_2_A", "2_1_B"]


def test_nodes_kinds_keys_and_labels():
    mol = make_mol(
        [
            protein_chain("A", 6),
            _lig(),
            _single("I", "ZN", "ZN", (3.8, 2, 2)),
        ]
    )
    asm = build_assembly(mol, EXCLUDED, min_heavy=6)
    assert asm.kinds == ["protein", "ligand", "ion"]
    assert asm.keys == ["P:AAAAAA", "L:XYZ", "I:ZN"]
    assert asm.labels == ["A", "L:XYZ", "I:ZN"]
    assert len(asm.xyz) == len(asm.node_idx) == len(asm.atom_ids) == 6 * 5 + 7 + 1
    assert asm.stats["n_merged"] == 0


def test_water_is_never_a_node_and_only_oxygens_are_returned():
    w = chain(
        "W",
        [residue("HOH", "X", [("O", "O", (5, 5, 5)), ("H1", "H", (5, 5, 6))])],
        etype="water",
    )
    asm = build_assembly(make_mol([protein_chain("A", 4), w]), EXCLUDED)
    assert asm.kinds == ["protein"]
    assert asm.stats["n_water"] == 1
    assert np.allclose(asm.water_xyz, [[5, 5, 5]])


def test_excluded_and_small_ligands_are_dropped_and_counted():
    mol = make_mol(
        [protein_chain("A", 4), _lig("E", "HEM_X"), _lig("S", "TINY", n=3), _lig("K")]
    )
    asm = build_assembly(mol, EXCLUDED, min_heavy=6)
    assert asm.keys[1:] == ["L:XYZ"]
    assert asm.stats["dropped"] == {"excluded": 1, "small": 1, "skipped_polymer": 0}
    assert asm.stats["excluded_names"] == {"HEM_X": 1}


def test_unknown_polymer_types_are_skipped():
    other = chain("P", protein_chain("P", 3)["residues"], etype="polysaccharide(D)")
    asm = build_assembly(make_mol([protein_chain("A", 4), other]), EXCLUDED)
    assert asm.kinds == ["protein"]
    assert asm.stats["dropped"]["skipped_polymer"] == 1


def test_nucleic_acid_chain_is_na_node_with_sequence_key():
    dna = chain(
        "D",
        [
            residue(c, o, [("C1'", "C", (3.0 * i, 0, 0)), ("N1", "N", (3.0 * i, 1, 0))])
            for i, (c, o) in enumerate([("DA", "A"), ("DC", "C"), ("DG", "G")])
        ],
        etype="polydeoxyribonucleotide",
    )
    asm = build_assembly(make_mol([dna]), EXCLUDED)
    assert asm.kinds == ["na"]
    assert asm.keys == ["N:ACG"]
    assert asm.donor.tolist() == [False, True] * 3


def test_branched_chain_is_a_ligand_named_after_its_sorted_components():
    atoms = lambda x: [(f"C{i}", "C", (x + i, 6.0, 6.0)) for i in range(4)]  # noqa: E731
    glycan = chain(
        "G",
        [residue("NAG", "X", atoms(0)), residue("BMA", "X", atoms(10))],
        etype="branched",
    )
    asm = build_assembly(make_mol([protein_chain("A", 3), glycan]), EXCLUDED)
    assert asm.kinds == ["protein", "ligand"]
    assert asm.keys[1] == "L:BMA-NAG"
    tiny = chain("T", [residue("NAG", "X", atoms(0)[:2])], etype="branched")
    asm = build_assembly(make_mol([protein_chain("A", 3), tiny]), EXCLUDED)
    assert asm.kinds == ["protein"]
    assert asm.stats["dropped"]["small"] == 1


def _atom(mol, chain_id_index):
    return int(mol.index_table.chains_to_atoms(np.array([chain_id_index]))[0])


def test_covalently_bonded_ligand_merges_into_the_polymer_node():
    n_prot_atoms = 4 * 5
    mol = make_mol(
        [protein_chain("A", 4), _lig()],
        struct_conn=[(1, n_prot_atoms, "covale")],
    )
    asm = build_assembly(mol, EXCLUDED)
    assert asm.kinds == ["protein"]
    assert asm.stats["n_merged"] == 1
    assert len(asm.xyz) == n_prot_atoms + 7


def test_metal_coordination_and_disulfide_never_merge_and_ions_never_merge():
    n_prot = 4 * 5
    mol = make_mol(
        [protein_chain("A", 4), _lig(), _single("I", "ZN", "ZN", (3.8, 2, 2))],
        struct_conn=[
            (1, n_prot, "metalc"),
            (1, n_prot, "disulf"),
            (1, n_prot + 7, "covale"),  # to an ion: never merged
        ],
    )
    asm = build_assembly(mol, EXCLUDED)
    assert asm.kinds == ["protein", "ligand", "ion"]
    assert asm.stats["n_merged"] == 0


def test_covalent_bond_between_two_polymers_does_not_merge_them():
    mol = make_mol(
        [protein_chain("A", 3), protein_chain("B", 3, z=4.0)],
        struct_conn=[(1, 15, "covale")],
    )
    assert build_assembly(mol, EXCLUDED).kinds == ["protein", "protein"]


def test_two_ligands_bonded_to_each_other_merge_to_the_lower_index():
    n_prot = 3 * 5
    mol = make_mol(
        [protein_chain("A", 3), _lig("L"), _lig("M", z=8.0)],
        struct_conn=[(n_prot + 7, n_prot, "covale")],
    )
    asm = build_assembly(mol, EXCLUDED)
    assert asm.kinds == ["protein", "ligand"]
    assert asm.stats["n_merged"] == 1
    assert (asm.node_idx == 1).sum() == 14


def test_entry_with_only_excluded_ligands_has_no_nodes():
    asm = build_assembly(make_mol([_lig("E", "HEM_X")]), EXCLUDED)
    assert asm.kinds == []
    assert len(asm.xyz) == 0
