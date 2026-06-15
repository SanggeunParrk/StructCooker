from pathlib import Path

import numpy as np
from datacooker import RecipeBook

from structcooker.instructions.transforms.openfold_structure import (
    assemble_cifmol,
    build_hierarchy,
    derive_atom_features,
    derive_bond_edges,
    derive_chain_features,
    derive_residue_features,
    load_ccd_entries,
)

"""Build an OpenFold3 distillation structure (CIFMol) Cooker.

Rebuilds the residue/chain hierarchy from the flat ``structure.npz`` atom table
and derives the chemistry fields from the CCD (keyed by residue name) so the
records match the existing cif LMDB content. Each conceptual step is its own
instruction.
"""

structure_recipe = RecipeBook()

structure_recipe.add(
    targets=(
        ("atom_to_res", np.ndarray),
        ("res_to_chain", np.ndarray),
        ("n_chain", int),
        ("res_names", np.ndarray),
        ("res_ids", np.ndarray),
        ("res_hetero", np.ndarray),
        ("chain_ids", np.ndarray),
        ("entity_ids", np.ndarray),
        ("chain_mol_types", np.ndarray),
    ),
    instruction=build_hierarchy,
    inputs={"kwargs": {"atom_arrays": ("atom_arrays", dict)}},
)

structure_recipe.add(
    targets=(("ccd_cache", dict),),
    instruction=load_ccd_entries,
    inputs={
        "kwargs": {
            "atom_arrays": ("atom_arrays", dict),
            "ccd_db_path": ("ccd_db_path", Path),
        },
    },
)

structure_recipe.add(
    targets=(("atom_features", dict),),
    instruction=derive_atom_features,
    inputs={"kwargs": {"atom_arrays": ("atom_arrays", dict), "ccd_cache": ("ccd_cache", dict)}},
)

structure_recipe.add(
    targets=(("bonds", dict),),
    instruction=derive_bond_edges,
    inputs={"kwargs": {"atom_arrays": ("atom_arrays", dict), "ccd_cache": ("ccd_cache", dict)}},
)

structure_recipe.add(
    targets=(("residue_features", dict),),
    instruction=derive_residue_features,
    inputs={
        "kwargs": {
            "res_names": ("res_names", np.ndarray),
            "res_ids": ("res_ids", np.ndarray),
            "res_hetero": ("res_hetero", np.ndarray),
            "ccd_cache": ("ccd_cache", dict),
        },
    },
)

structure_recipe.add(
    targets=(("chain_features", dict),),
    instruction=derive_chain_features,
    inputs={
        "kwargs": {
            "chain_ids": ("chain_ids", np.ndarray),
            "entity_ids": ("entity_ids", np.ndarray),
            "chain_mol_types": ("chain_mol_types", np.ndarray),
        },
    },
)

structure_recipe.add(
    targets=(("cifmol_dict", dict),),
    instruction=assemble_cifmol,
    inputs={
        "kwargs": {
            "atom_arrays": ("atom_arrays", dict),
            "atom_to_res": ("atom_to_res", np.ndarray),
            "res_to_chain": ("res_to_chain", np.ndarray),
            "n_chain": ("n_chain", int),
            "atom_features": ("atom_features", dict),
            "bonds": ("bonds", dict),
            "residue_features": ("residue_features", dict),
            "chain_features": ("chain_features", dict),
        },
    },
)

RECIPE = structure_recipe
TARGETS = ["cifmol_dict"]
