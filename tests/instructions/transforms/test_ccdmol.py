import numpy as np
import pytest

from structcooker.instructions.transforms.ccdmol import convert_chem_comp_to_ccdmol


def _record(n_atoms, bonds):
    src, dst = (np.array([b[0] for b in bonds]), np.array([b[1] for b in bonds])) if bonds else (np.array([]), np.array([]))
    nodes = {k: {"value": np.array(v)} for k, v in {
        "id": [f"A{i}" for i in range(n_atoms)], "element": ["C"] * n_atoms, "aromatic": [False] * n_atoms,
        "stereo": ["N"] * n_atoms, "charge": [0] * n_atoms, "model_xyz": np.zeros((n_atoms, 3))}.items()}
    edges = {k: {"value": np.array([1] * len(bonds)), "src_indices": src, "dst_indices": dst}
             for k in ("bond_type", "bond_aromatic", "bond_stereo")}
    residue = {"nodes": {"id": {"value": np.array(["X"])}, "formula": {"value": np.array(["C2"])}}}
    return {"chem_comp_dict": {"atom": {"nodes": nodes, "edges": edges}, "residue": residue}}


def test_reshape_keeps_primary_fields_and_one_residue_chain():
    mol = convert_chem_comp_to_ccdmol(_record(2, [(0, 1)]), "ETH")
    assert list(mol["atoms"]["nodes"]["id"]["value"]) == ["A0", "A1"]
    assert list(mol["atoms"]["edges"]["bond_type"]["src_indices"]) == [0]
    assert list(mol["chains"]["nodes"]["id"]["value"]) == ["ETH"]
    assert list(mol["index_table"]["res_atom_indptr"]) == [0, 2]


def test_single_atom_component_passes_and_atomless_is_refused():
    ion = convert_chem_comp_to_ccdmol(_record(1, []), "ZN")
    assert len(ion["atoms"]["edges"]["bond_type"]["value"]) == 0
    with pytest.raises(ValueError, match="no atoms"):
        convert_chem_comp_to_ccdmol({"chem_comp_dict": None}, "UNL")
