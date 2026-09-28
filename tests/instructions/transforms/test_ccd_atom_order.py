import numpy as np
from biomol.core.container import FeatureContainer
from biomol.core.feature import EdgeFeature, NodeFeature

from structcooker.instructions.transforms.cif import _compare_each_chem_comp


def test_inherited_ccd_features_follow_atom_ids_and_remap_bonds():
    residue = FeatureContainer({"id": NodeFeature(np.array(["LIG"]))})
    cif = {"residue": residue, "atom": FeatureContainer({"id": NodeFeature(np.array(["B", "A"]))})}
    ideal = {"residue": residue, "atom": FeatureContainer({
        "id": NodeFeature(np.array(["A", "B"])),
        "model_xyz": NodeFeature(np.array([[1., 2., 3.], [4., 5., 6.]])),
        "bond": EdgeFeature(np.array(["SING"]), src_indices=np.array([0]), dst_indices=np.array([1])),
    })}
    merged = _compare_each_chem_comp(cif, ideal)["atom"]
    np.testing.assert_array_equal(merged["model_xyz"].value, [[4., 5., 6.], [1., 2., 3.]])
    np.testing.assert_array_equal(merged["bond"].src_indices, [1])
    np.testing.assert_array_equal(merged["bond"].dst_indices, [0])
    np.testing.assert_array_equal(ideal["atom"]["id"].value, ["A", "B"])


def test_missing_cif_atom_requires_matching_reference(monkeypatch):
    from pathlib import Path

    import pytest

    from structcooker.instructions.transforms.cif import compare_chem_comp

    residue = FeatureContainer({"id": NodeFeature(np.array(["LIG"]))})
    cif = {"LIG": {"residue": residue, "atom": FeatureContainer({
        "id": NodeFeature(np.array(["A", "MG1"])),
    })}}
    primary = {"chem_comp_dict": {"residue": residue.to_dict(), "atom": FeatureContainer({
        "id": NodeFeature(np.array(["A"])),
        "charge": NodeFeature(np.array(["0"])),
    }).to_dict()}}
    monkeypatch.setattr("structcooker.instructions.transforms.cif.read_lmdb", lambda *_: primary)
    with pytest.raises(ValueError, match="lacks CIF atoms.*MG1"):
        compare_chem_comp(cif, Path("primary"))
    with pytest.raises(ValueError, match="lacks CIF atoms.*MG1"):
        compare_chem_comp(cif, Path("primary"), Path("nonmatching"))


def test_reference_preserves_cif_atoms_and_aligns_inherited_chemistry(monkeypatch):
    from pathlib import Path

    from structcooker.instructions.transforms.cif import compare_chem_comp

    residue = FeatureContainer({"id": NodeFeature(np.array(["LIG"]))})
    cif = {"LIG": {"residue": residue, "atom": FeatureContainer({
        "id": NodeFeature(np.array(["A", "MG1"])),
    })}}
    def component(ids, charges):
        return {"chem_comp_dict": {"residue": residue.to_dict(), "atom": FeatureContainer({
            "id": NodeFeature(np.array(ids)), "charge": NodeFeature(np.array(charges)),
        }).to_dict()}}
    sources = {Path("primary"): component(["A"], ["0"]),
               Path("reference"): component(["MG1", "A"], ["2", "0"])}
    monkeypatch.setattr("structcooker.instructions.transforms.cif.read_lmdb", lambda path, _: sources[path])
    result = compare_chem_comp(cif, Path("primary"), Path("reference"))["LIG"]["atom"]
    np.testing.assert_array_equal(result["id"].value, ["A", "MG1"])
    np.testing.assert_array_equal(result["charge"].value, ["0", "2"])


def test_exact_ccd_aliases_preserve_cif_names_without_reference(monkeypatch):
    from pathlib import Path

    from structcooker.instructions.transforms.cif import compare_chem_comp

    residue = FeatureContainer({"id": NodeFeature(np.array(["LIG"]))})
    cif = {"LIG": {"residue": residue, "atom": FeatureContainer({
        "id": NodeFeature(np.array(["C1", "N1"])),
    })}}
    primary = {"chem_comp_dict": {"residue": residue.to_dict(), "atom": FeatureContainer({
        "id": NodeFeature(np.array(["NAA", "CAB"])),
        "alt_atom_id": NodeFeature(np.array(["N1", "C1"])),
        "charge": NodeFeature(np.array(["1", "0"])),
    }).to_dict()}}
    monkeypatch.setattr("structcooker.instructions.transforms.cif.read_lmdb", lambda *_: primary)
    result = compare_chem_comp(cif, Path("primary"))["LIG"]["atom"]
    np.testing.assert_array_equal(result["id"].value, ["C1", "N1"])
    np.testing.assert_array_equal(result["charge"].value, ["0", "1"])
    assert "alt_atom_id" not in result


def test_component_atom_names_are_supported_and_ambiguous_aliases_rejected():
    import pytest

    from structcooker.instructions.transforms.cif import _match_ccd_atom_aliases

    nodes = {"id": {"value": np.array(["CAA", "NBB"])},
             "component_atom_id": {"value": np.array(["C1", "N1"])}}
    _match_ccd_atom_aliases({"atom": {"nodes": nodes}}, {"C1", "N1"})
    np.testing.assert_array_equal(nodes["id"]["value"], ["C1", "N1"])
    nodes["id"] = {"value": np.array(["CAA", "NBB"])}
    nodes["alt_atom_id"] = {"value": np.array(["N1", "C1"])}
    with pytest.raises(ValueError, match="alias tables disagree"):
        _match_ccd_atom_aliases({"atom": {"nodes": nodes}}, {"C1", "N1"})
