"""Rebuild a (rewrapped) attached CIF LMDB to add the chain contact graph.

The disordered CIF was produced by ``rewrap_cif`` (pure re-layout of the
OpenFold records), so it never went through the full ``cif.py`` ingest that runs
``extract_contact_graph``. As a result ``chains.contact`` is empty and every
entry looks monomeric to ``extract_edge_node``. This recipe recomputes the
6 A chain-chain contact graph from the coordinates already present in each
attached CIFMol and writes it back under the same ``cifmol_attached_dict``
layout, so a subsequent ``edge_node`` extract recovers the interfaces.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from biomol.core.feature import EdgeFeature
from biomol.enums import StructureLevel
from datacooker import RecipeBook

from structcooker.instructions.transforms.geometry import chain_contacts_grid
from structcooker.mols import CIFMolAttached

if TYPE_CHECKING:
    from biomol.core.types import BioMolDict

_D_THR = 6.0


def attach_chain_contacts(cifmol: CIFMolAttached) -> BioMolDict:
    """Recompute the chain contact graph for one attached CIFMol.

    Mirrors ``geometry.extract_contact_graph`` (6 A threshold) but operates on an
    already-built ``CIFMolAttached`` and returns the storage dict so the rebuild
    round-trips into the same ``cifmol_attached_dict`` layout.
    """
    xyz = cifmol.atoms.xyz.value
    chain_idx = cifmol.index_table.atoms_to_chains(np.arange(xyz.shape[0]))
    src, dst, counts = chain_contacts_grid(xyz, chain_idx, _D_THR, count_atom_pairs_once=True)
    contact = EdgeFeature(
        value=counts.astype(np.int32),
        src_indices=src.astype(chain_idx.dtype),
        dst_indices=dst.astype(chain_idx.dtype),
    )
    updated = cifmol.update_features(StructureLevel.CHAIN, contact=contact)
    return updated.to_dict()


recipe = RecipeBook()

recipe.step(
    outputs=(("cifmol_attached_dict", dict),),
    instruction=attach_chain_contacts,
    kwargs={"cifmol": ("cifmol", CIFMolAttached)},
)

RECIPE = recipe
TARGETS = ["cifmol_attached_dict"]
