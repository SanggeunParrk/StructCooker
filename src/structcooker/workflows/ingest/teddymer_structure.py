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
    wrap_cifmol,
)
from structcooker.instructions.transforms.teddymer import attach_bonds

"""Build a Teddymer dimer structure (CIFMol) Cooker.

Identical to ``openfold_structure`` from ``atom_site_dict`` onward -- the same hierarchy
rebuild and the same CCD-derived chemistry -- with one step in front. A PDB carries no
connectivity, so ``attach_bonds`` derives the atom-index pairs from the CCD and produces
the ``atom_site_dict`` those steps expect out of the reader's ``raw_atom_site_dict``.

``ccd_cache`` is keyed by residue name only, so it is loaded from the raw table; that is
what keeps the bond step from depending on its own output.
"""

teddymer_recipe = RecipeBook()

teddymer_recipe.add(
    targets=(("ccd_cache", dict),),
    instruction=load_ccd_entries,
    inputs={
        "kwargs": {
            "atom_site_dict": ("raw_atom_site_dict", dict),
            "ccd_db_path": ("ccd_db_path", Path),
        },
    },
)

teddymer_recipe.add(
    targets=(("atom_site_dict", dict),),
    instruction=attach_bonds,
    inputs={
        "kwargs": {
            "raw_atom_site_dict": ("raw_atom_site_dict", dict),
            "ccd_cache": ("ccd_cache", dict),
        },
    },
)

teddymer_recipe.add(
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
    inputs={"kwargs": {"atom_site_dict": ("atom_site_dict", dict)}},
)

teddymer_recipe.add(
    targets=(("atom_features", dict),),
    instruction=derive_atom_features,
    inputs={"kwargs": {"atom_site_dict": ("atom_site_dict", dict), "ccd_cache": ("ccd_cache", dict)}},
)

teddymer_recipe.add(
    targets=(("bonds", dict),),
    instruction=derive_bond_edges,
    inputs={"kwargs": {"atom_site_dict": ("atom_site_dict", dict), "ccd_cache": ("ccd_cache", dict)}},
)

teddymer_recipe.add(
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

teddymer_recipe.add(
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

teddymer_recipe.add(
    targets=(("cifmol_dict", dict),),
    instruction=assemble_cifmol,
    inputs={
        "kwargs": {
            "atom_site_dict": ("atom_site_dict", dict),
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

teddymer_recipe.add(
    targets=(("assembly_dict", dict), ("metadata_dict", dict)),
    instruction=wrap_cifmol,
    inputs={
        "kwargs": {
            "cifmol_dict": ("cifmol_dict", dict),
            "entry_id": ("entry_id", str),
        },
    },
)

RECIPE = teddymer_recipe
TARGETS = ["assembly_dict", "metadata_dict"]
