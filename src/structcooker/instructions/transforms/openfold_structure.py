"""Build unified ``CIFMol`` payloads from OpenFold3 distillation structures.

``structure.npz`` only stores a flat per-atom table, so the residue/chain
hierarchy is reconstructed here and every feature absent from the distillation
output (``aromatic``, ``stereo``, ``b_factor``, ``formula``,
``one_letter_code``, ``entity_type`` ...) is materialised as an empty field so
the records share the schema of the existing ``cif.lmdb`` BioMol database.
"""

from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np
from biomol.core import EdgeFeature, FeatureContainer, NodeFeature
from biomol.core.index import IndexTable

from structcooker.mols import CIFMol

if TYPE_CHECKING:
    from biomol.core.types import BioMolDict


def _group_ids(*columns: np.ndarray) -> np.ndarray:
    """Assign a contiguous group id whenever any column value changes.

    Atoms in ``structure.npz`` are laid out residue-by-residue (and residues
    chain-by-chain), so contiguous runs of equal keys form one parent.
    """
    n = len(columns[0])
    group_ids = np.zeros(n, dtype=np.int64)
    for i in range(1, n):
        changed = any(col[i] != col[i - 1] for col in columns)
        group_ids[i] = group_ids[i - 1] + (1 if changed else 0)
    return group_ids


def _first_per_group(group_ids: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Return the value of the first member of each contiguous group."""
    _, first_idx = np.unique(group_ids, return_index=True)
    return values[first_idx]


def _empty_str(shape: tuple[int, ...]) -> np.ndarray:
    """Return a placeholder string array for a missing feature."""
    return np.full(shape, "", dtype="<U1")


def _empty_edge(value_dtype: str = "<U1") -> EdgeFeature:
    """Return an edge feature with no edges (a materialised but empty field)."""
    return EdgeFeature(
        value=np.empty((0,), dtype=value_dtype),
        src_indices=np.empty((0,), dtype=np.int64),
        dst_indices=np.empty((0,), dtype=np.int64),
    )


def _build_cifmol(atoms: dict[str, np.ndarray]) -> CIFMol:
    """Reconstruct a ``CIFMol`` from a flat per-atom feature table."""
    n_atoms = len(atoms["coord"])
    chain_id = atoms["chain_id"].astype(str)
    res_id = atoms["res_id"]
    ins_code = atoms["ins_code"].astype(str)
    entity_id = atoms["entity_id"].astype(str)

    atom_to_res = _group_ids(chain_id, res_id, ins_code)
    res_chain_id = _first_per_group(atom_to_res, chain_id)
    res_entity_id = _first_per_group(atom_to_res, entity_id)
    res_to_chain = _group_ids(res_chain_id, res_entity_id)

    bonds = atoms["bonds"]
    bond_src = bonds[:, 0].astype(np.int64)
    bond_dst = bonds[:, 1].astype(np.int64)
    bond_order = bonds[:, 2].astype(str).astype("<U21")
    n_bonds = len(bonds)

    atom_container = FeatureContainer(
        {
            "id": NodeFeature(atoms["atom_name"].astype("<U4")),
            "element": NodeFeature(atoms["element"].astype("<U2")),
            "aromatic": NodeFeature(_empty_str((n_atoms,))),
            "stereo": NodeFeature(_empty_str((n_atoms,))),
            "charge": NodeFeature(atoms["charge"].astype(str)),
            "model_xyz": NodeFeature(np.full((n_atoms, 3), "", dtype="<U7")),
            "xyz": NodeFeature(atoms["coord"].astype(np.float64)),
            "b_factor": NodeFeature(np.zeros(n_atoms, dtype=np.float64)),
            "occupancy": NodeFeature(atoms["occupancy"].astype(np.float64)),
            "bond_type": EdgeFeature(bond_order, bond_src, bond_dst),
            "bond_aromatic": EdgeFeature(_empty_str((n_bonds,)), bond_src, bond_dst),
            "bond_stereo": EdgeFeature(_empty_str((n_bonds,)), bond_src, bond_dst),
            "struct_conn": _empty_edge("<U6"),
        },
    )

    res_name = _first_per_group(atom_to_res, atoms["res_name"].astype(str))
    res_res_id = _first_per_group(atom_to_res, res_id).astype(str)
    res_hetero = _first_per_group(atom_to_res, atoms["hetero"]).astype(np.int64)
    n_res = len(res_name)
    residue_container = FeatureContainer(
        {
            "id": NodeFeature(res_res_id.astype("<U34")),
            "formula": NodeFeature(_empty_str((n_res,))),
            "cif_idx": NodeFeature(res_res_id.astype("<U21")),
            "auth_idx": NodeFeature(res_res_id.astype("<U2")),
            "chem_comp_id": NodeFeature(res_name.astype("<U3")),
            "hetero": NodeFeature(res_hetero),
            "one_letter_code_can": NodeFeature(_empty_str((n_res,))),
            "one_letter_code": NodeFeature(_empty_str((n_res,))),
            "bond": _empty_edge("<U21"),
            "struct_conn": _empty_edge("int64"),
        },
    )

    chain_chain_id = _first_per_group(res_to_chain, res_chain_id)
    chain_entity_id = _first_per_group(res_to_chain, res_entity_id)
    n_chain = len(chain_chain_id)
    chain_container = FeatureContainer(
        {
            "entity_id": NodeFeature(chain_entity_id.astype("<U1")),
            "entity_type": NodeFeature(_empty_str((n_chain,))),
            "auth_asym_id": NodeFeature(_empty_str((n_chain,))),
            "chain_id": NodeFeature(chain_chain_id.astype("<U3")),
            "contact": _empty_edge("int32"),
        },
    )

    index_table = IndexTable.from_parents(atom_to_res, res_to_chain, n_chain=n_chain)
    return CIFMol(
        atom_container=atom_container,
        residue_container=residue_container,
        chain_container=chain_container,
        index_table=index_table,
    )


def build_structure_dict() -> Callable[..., "BioMolDict"]:
    """Reconstruct a unified ``CIFMol`` payload from the flat atom arrays."""

    def _worker(atom_arrays: dict[str, np.ndarray]) -> "BioMolDict":
        return _build_cifmol(atom_arrays).to_dict()

    return _worker
