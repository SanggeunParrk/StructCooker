import numpy as np
import pytest

from structcooker.instructions.transforms.cif import _rearrange_each_asym_id


@pytest.mark.parametrize("size", [0, 1, 64, 200])
@pytest.mark.parametrize("remove_hydrogen", [False, True])
def test_altloc_grouping_matches_original_per_atom_selection(size, remove_hydrogen):
    rng = np.random.default_rng(20260910 + size)
    def strings(values):
        return np.asarray(values, dtype=str)
    data = {
        "label_seq_id": strings(rng.integers(1, 5, size)),
        "pdbx_PDB_ins_code": rng.choice([".", "?", "A"], size),
        "auth_seq_id": strings(rng.integers(1, 5, size)),
        "Cartn_x": strings(np.arange(size)),
        "Cartn_y": strings(np.arange(size) + 1),
        "Cartn_z": strings(np.arange(size) + 2),
        "B_iso_or_equiv": strings(np.ones(size)),
        "occupancy": strings(np.ones(size)),
        "type_symbol": rng.choice(["C", "N", "H", "D"], size),
        "label_atom_id": rng.choice(["CA", "N", "C"], size),
        "label_comp_id": np.full(size, "ALA"),
        "label_alt_id": rng.choice([".", "A", "B"], size),
        "pdbx_PDB_model_num": rng.choice(["1", "2", "10"], size),
        "auth_asym_id": np.full(size, "A"),
    }
    keep = ~np.isin(data["type_symbol"], ["H", "D"]) if remove_hydrogen else np.ones(size, dtype=bool)
    raw = {k: v[keep] for k, v in data.items()}
    auth = [a if ins in {".", "?"} else f"{a}.{ins}"
            for a, ins in zip(raw["auth_seq_id"], raw["pdbx_PDB_ins_code"], strict=True)]
    determinants = np.asarray([f"{atom}.{res}" for atom, res in zip(raw["label_atom_id"], auth, strict=True)])
    result = _rearrange_each_asym_id(dict(data), remove_hydrogen)["atom_site"]
    assert set(result) == set(raw["pdbx_PDB_model_num"])
    for model in np.unique(raw["pdbx_PDB_model_num"]):
        model_mask = raw["pdbx_PDB_model_num"] == model
        for alt in np.unique(raw["label_alt_id"][model_mask]):
            selected = np.zeros(len(determinants), dtype=bool)
            for atom in np.unique(determinants):
                group = determinants == atom
                target = alt if np.any(raw["label_alt_id"][group] == alt) else "."
                selected |= group & (raw["label_alt_id"] == target)
            expected = selected & model_mask
            actual = result[model][alt]
            np.testing.assert_array_equal(actual["xyz"][:, 0], raw["Cartn_x"][expected].astype(float))
            np.testing.assert_array_equal(actual["atom"], raw["label_atom_id"][expected])
            np.testing.assert_array_equal(actual["auth_idx"], np.asarray(auth)[expected])
