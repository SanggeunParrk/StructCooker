"""Assemble a unified ``CIFMol`` from OpenFold3 distillation structures.

``structure.npz`` only ships a flat per-atom table, so the residue/chain
hierarchy is rebuilt and the chemistry fields the distillation output omits
(``aromatic`` / ``stereo`` / ``charge`` / ``model_xyz`` / ``formula`` /
bond types) are **derived from the CCD** keyed by residue name, exactly as the
mmCIF ingest does. ``entity_type`` / ``one_letter_code`` are derived from
``molecule_type_id`` / a residue map. The work is split into one instruction
per conceptual step (group, CCD lookup, per-level feature derivation, assemble).
"""

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from biomol.core import EdgeFeature, FeatureContainer, NodeFeature
from biomol.core.index import IndexTable

from structcooker.instructions.readers.lmdb import read_lmdb
from structcooker.mols import CIFMol
from structcooker.utils.mapping import MoleculeType

if TYPE_CHECKING:
    from biomol.core.types import BioMolDict

_PROTEIN_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
_RNA_1 = {"A": "A", "C": "C", "G": "G", "U": "U"}
# OpenFold3 distillation molecule_type_id encoding (protein=0, rna=1, ...),
# which differs from StructCooker's internal MoleculeType ordering.
_OPENFOLD_MOL_TYPE = {
    0: MoleculeType.PROTEIN,
    1: MoleculeType.RNA,
    2: MoleculeType.DNA,
    3: MoleculeType.LIGAND,
}


def _group_ids(*columns: np.ndarray) -> np.ndarray:
    """Assign a contiguous group id whenever any column value changes."""
    n = len(columns[0])
    group_ids = np.zeros(n, dtype=np.int64)
    for i in range(1, n):
        if any(col[i] != col[i - 1] for col in columns):
            group_ids[i] = group_ids[i - 1] + 1
        else:
            group_ids[i] = group_ids[i - 1]
    return group_ids


def _first_per_group(group_ids: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Return the value of the first member of each contiguous group."""
    _, first_idx = np.unique(group_ids, return_index=True)
    return values[first_idx]


def _empty_edge(value_shape: tuple[int, ...], value_dtype: str) -> EdgeFeature:
    """Return an edge feature with no edges (a materialised but empty field)."""
    return EdgeFeature(
        value=np.empty(value_shape, dtype=value_dtype),
        src_indices=np.empty((0,), dtype=np.int64),
        dst_indices=np.empty((0,), dtype=np.int64),
    )


def _read_ccd_entry(ccd_db_path: Path, res_name: str) -> dict[str, Any] | None:
    """Read one CCD component, returning ``None`` if it is missing."""
    try:
        return read_lmdb(ccd_db_path, res_name)["chem_comp_dict"]
    except Exception:  # noqa: BLE001 - missing component → derive defaults
        return None


def _ccd_atom_lookup(entry: dict[str, Any] | None) -> dict[str, dict[str, Any]]:
    """Map atom name → CCD atom features for one component."""
    if entry is None:
        return {}
    nodes = entry["atom"]["nodes"]
    names = nodes["id"]["value"]
    return {
        name: {
            "aromatic": str(nodes["aromatic"]["value"][i]),
            "stereo": str(nodes["stereo"]["value"][i]),
            "charge": str(nodes["charge"]["value"][i]),
            "model_xyz": nodes["model_xyz"]["value"][i],
        }
        for i, name in enumerate(names.tolist())
    }


def build_hierarchy(atom_arrays: dict[str, np.ndarray]) -> dict[str, Any]:
    """Recover the atom→residue→chain hierarchy from the flat atom table."""
    chain_id = atom_arrays["chain_id"].astype(str)
    res_id = atom_arrays["res_id"]
    ins_code = atom_arrays["ins_code"].astype(str)
    entity_id = atom_arrays["entity_id"].astype(str)
    mol_type = atom_arrays["molecule_type_id"]

    atom_to_res = _group_ids(chain_id, res_id, ins_code)
    res_chain_id = _first_per_group(atom_to_res, chain_id)
    res_entity_id = _first_per_group(atom_to_res, entity_id)
    res_to_chain = _group_ids(res_chain_id, res_entity_id)
    return {
        "atom_to_res": atom_to_res,
        "res_to_chain": res_to_chain,
        "n_chain": int(res_to_chain[-1] + 1) if len(res_to_chain) else 0,
        "res_name": _first_per_group(atom_to_res, atom_arrays["res_name"].astype(str)),
        "res_id": _first_per_group(atom_to_res, res_id),
        "res_hetero": _first_per_group(atom_to_res, atom_arrays["hetero"]).astype(np.int64),
        "chain_id": _first_per_group(res_to_chain, res_chain_id),
        "entity_id": _first_per_group(res_to_chain, res_entity_id),
        "chain_mol_type": _first_per_group(
            res_to_chain, _first_per_group(atom_to_res, mol_type),
        ),
    }


def load_ccd_entries(
    atom_arrays: dict[str, np.ndarray],
    ccd_db_path: Path,
) -> dict[str, Any]:
    """Load the CCD component for every residue name present in the structure."""
    path = Path(ccd_db_path)
    names = np.unique(atom_arrays["res_name"].astype(str)).tolist()
    return {res_name: _read_ccd_entry(path, res_name) for res_name in names}


def derive_atom_features(
    atom_arrays: dict[str, np.ndarray],
    ccd_cache: dict[str, Any],
) -> dict[str, np.ndarray]:
    """Fill per-atom ``aromatic`` / ``stereo`` / ``charge`` / ``model_xyz`` from the CCD."""
    res_name = atom_arrays["res_name"].astype(str)
    atom_name = atom_arrays["atom_name"].astype(str)
    n = len(atom_name)
    aromatic = np.full(n, "N", dtype="<U1")
    stereo = np.full(n, "N", dtype="<U1")
    charge = np.full(n, "0", dtype="<U2")
    model_xyz = np.full((n, 3), "", dtype="<U7")
    lookups = {r: _ccd_atom_lookup(ccd_cache.get(r)) for r in np.unique(res_name).tolist()}
    for i in range(n):
        feat = lookups[res_name[i]].get(atom_name[i])
        if feat is None:
            continue
        aromatic[i] = feat["aromatic"]
        stereo[i] = feat["stereo"]
        charge[i] = feat["charge"]
        model_xyz[i] = feat["model_xyz"]
    return {"aromatic": aromatic, "stereo": stereo, "charge": charge, "model_xyz": model_xyz}


def derive_bond_edges(
    atom_arrays: dict[str, np.ndarray],
    ccd_cache: dict[str, Any],
) -> dict[str, np.ndarray]:
    """Build atom- and residue-level bond edges from connectivity + CCD chemistry."""
    res_name = atom_arrays["res_name"].astype(str)
    atom_name = atom_arrays["atom_name"].astype(str)
    atom_to_res = _group_ids(
        atom_arrays["chain_id"].astype(str),
        atom_arrays["res_id"],
        atom_arrays["ins_code"].astype(str),
    )
    # CCD intra-residue bond chemistry keyed by (res_name, frozenset{atom_a, atom_b})
    chem: dict[tuple[str, frozenset[str]], tuple[str, str, str]] = {}
    for r, entry in ccd_cache.items():
        if entry is None:
            continue
        edges = entry["atom"]["edges"]
        names = entry["atom"]["nodes"]["id"]["value"]
        bt = edges["bond_type"]
        for k in range(len(bt["value"])):
            a = str(names[bt["src_indices"][k]])
            b = str(names[bt["dst_indices"][k]])
            chem[(r, frozenset({a, b}))] = (
                str(bt["value"][k]),
                str(edges["bond_aromatic"]["value"][k]),
                str(edges["bond_stereo"]["value"][k]),
            )
    bonds = atom_arrays["bonds"]
    src = bonds[:, 0].astype(np.int64)
    dst = bonds[:, 1].astype(np.int64)
    m = len(bonds)
    bond_type = np.full(m, "SING", dtype="<U21")
    bond_aromatic = np.full(m, "N", dtype="<U1")
    bond_stereo = np.full(m, "N", dtype="<U1")
    res_bonds: dict[tuple[int, int], str] = {}
    for k in range(m):
        i, j = src[k], dst[k]
        if atom_to_res[i] == atom_to_res[j]:
            hit = chem.get(
                (str(res_name[i]), frozenset({str(atom_name[i]), str(atom_name[j])})),
            )
            if hit is not None:
                bond_type[k], bond_aromatic[k], bond_stereo[k] = hit
        else:
            res_bonds.setdefault((int(atom_to_res[i]), int(atom_to_res[j])), str(bond_type[k]))
    return {
        "src": src,
        "dst": dst,
        "bond_type": bond_type,
        "bond_aromatic": bond_aromatic,
        "bond_stereo": bond_stereo,
        "res_src": np.array([p[0] for p in res_bonds], dtype=np.int64),
        "res_dst": np.array([p[1] for p in res_bonds], dtype=np.int64),
        "res_bond_value": np.array(list(res_bonds.values()), dtype="<U21"),
    }


def derive_residue_features(
    hierarchy: dict[str, Any],
    ccd_cache: dict[str, Any],
) -> dict[str, np.ndarray]:
    """Derive residue-level features (formula from CCD, one-letter codes, indices)."""
    res_name = hierarchy["res_name"]
    res_id = hierarchy["res_id"].astype(str)
    n = len(res_name)
    formula = np.full(n, "", dtype="<U15")
    one_letter = np.full(n, "X", dtype="<U1")
    for i, r in enumerate(res_name.tolist()):
        entry = ccd_cache.get(r)
        if entry is not None:
            formula[i] = str(entry["residue"]["nodes"]["formula"]["value"][0])[:15]
        one_letter[i] = _PROTEIN_3TO1.get(r, _RNA_1.get(r, "X"))
    return {
        "id": res_name.astype("<U34"),
        "formula": formula,
        "cif_idx": res_id.astype("<U21"),
        "auth_idx": res_id.astype("<U21"),
        "chem_comp_id": res_name.astype("<U5"),
        "hetero": hierarchy["res_hetero"],
        "one_letter_code_can": one_letter,
        "one_letter_code": one_letter,
    }


def derive_chain_features(hierarchy: dict[str, Any]) -> dict[str, np.ndarray]:
    """Derive chain-level features (``entity_type`` from ``molecule_type_id``)."""
    chain_id = hierarchy["chain_id"].astype(str)
    entity_id = hierarchy["entity_id"].astype(str)
    entity_type = np.array(
        [
            _OPENFOLD_MOL_TYPE.get(int(t), MoleculeType.NA).value
            for t in hierarchy["chain_mol_type"]
        ],
        dtype="<U49",
    )
    return {
        "entity_id": entity_id.astype("<U8"),
        "entity_type": entity_type,
        "auth_asym_id": chain_id.astype("<U8"),
        "chain_id": chain_id.astype("<U8"),
    }


def assemble_cifmol(  # noqa: PLR0913
    atom_arrays: dict[str, np.ndarray],
    hierarchy: dict[str, Any],
    atom_features: dict[str, np.ndarray],
    bonds: dict[str, np.ndarray],
    residue_features: dict[str, np.ndarray],
    chain_features: dict[str, np.ndarray],
) -> "BioMolDict":
    """Assemble the derived features into a unified ``CIFMol`` payload."""
    n_atoms = len(atom_arrays["coord"])
    atoms = FeatureContainer(
        {
            "id": NodeFeature(atom_arrays["atom_name"].astype("<U4")),
            "element": NodeFeature(atom_arrays["element"].astype("<U2")),
            "aromatic": NodeFeature(atom_features["aromatic"]),
            "stereo": NodeFeature(atom_features["stereo"]),
            "charge": NodeFeature(atom_features["charge"]),
            "model_xyz": NodeFeature(atom_features["model_xyz"]),
            "xyz": NodeFeature(atom_arrays["coord"].astype(np.float64)),
            "b_factor": NodeFeature(np.zeros(n_atoms, dtype=np.float64)),
            "occupancy": NodeFeature(atom_arrays["occupancy"].astype(np.float64)),
            "bond_type": EdgeFeature(bonds["bond_type"], bonds["src"], bonds["dst"]),
            "bond_aromatic": EdgeFeature(bonds["bond_aromatic"], bonds["src"], bonds["dst"]),
            "bond_stereo": EdgeFeature(bonds["bond_stereo"], bonds["src"], bonds["dst"]),
            "struct_conn": _empty_edge((0, 2), "<U6"),
        },
    )
    residue_fields: dict[str, NodeFeature | EdgeFeature] = {
        k: NodeFeature(v) for k, v in residue_features.items()
    }
    residue_fields["bond"] = EdgeFeature(
        bonds["res_bond_value"], bonds["res_src"], bonds["res_dst"],
    )
    residue_fields["struct_conn"] = _empty_edge((0,), "int64")
    chain_fields: dict[str, NodeFeature | EdgeFeature] = {
        k: NodeFeature(v) for k, v in chain_features.items()
    }
    chain_fields["contact"] = _empty_edge((0,), "int32")
    index_table = IndexTable.from_parents(
        hierarchy["atom_to_res"], hierarchy["res_to_chain"], n_chain=hierarchy["n_chain"],
    )
    return CIFMol(
        atom_container=atoms,
        residue_container=FeatureContainer(residue_fields),
        chain_container=FeatureContainer(chain_fields),
        index_table=index_table,
    ).to_dict()
