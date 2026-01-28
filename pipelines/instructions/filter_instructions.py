from collections.abc import Callable
from datetime import date
from typing import cast

import numpy as np
from biomol.core import FeatureContainer, NodeFeature
from biomol.core.types import FeatureContainerDict

from pipelines.cifmol import CIFMol
from pipelines.instructions.seq_instructions import extract_sequence_from_cifmol


def filter_by_resolution_and_date(
    resolution_cutoff: float = 9.0,
    date_cutoff: date = date(2099, 1, 1),
) -> Callable[[CIFMol|None], CIFMol|None]:
    """Filter instruction to select entries by resolution and date."""
    def worker(cifmol: CIFMol|None) -> CIFMol|None:
        if cifmol is None:
            return None

        resolution, deposition_date = cifmol.metadata["resolution"], cifmol.metadata["deposition_date"]
        deposition_date = date.fromisoformat(deposition_date)
        if resolution is not None and resolution <= resolution_cutoff and deposition_date < date_cutoff:
            return cifmol
        return None

    return worker

def filter_water(cifmol: CIFMol|None) -> CIFMol|None:
    """Filter instruction to remove water molecules from CIFMol."""
    if cifmol is None:
        return None
    water_mask = ~np.isin(cifmol.residues.chem_comp_id, ["HOH", "DOD"])
    if water_mask.sum() == 0:
        return None
    cifmol = cifmol.residues[water_mask].extract()

    if len(cifmol.chains) == 0:
        return None

    return cifmol

def filter_signalp(
    cifmol: CIFMol|None,
    seq2seqid_dict:dict[str,str],
    signalp_dict:dict[str,tuple[int,int]],
) -> dict|None:
    """Filter instruction to remove signal peptides from CIFMol."""
    if cifmol is None:
        return None
    seq_dict = extract_sequence_from_cifmol(cifmol)
    valid_residue_indices = []
    cursor = 0
    for chain_id, seq in seq_dict.items():
        seq_id = seq2seqid_dict[seq]
        chain_cifmol = cifmol.chains[cifmol.chains.chain_id == chain_id].extract()
        if seq_id not in signalp_dict:
            valid_residue_indices.extend(list(range(cursor, cursor + len(chain_cifmol.residues))))
            cursor += len(chain_cifmol.residues)
            continue
        _, signalp_end = signalp_dict[seq_id]
        valid_residue_indices.extend(list(range(cursor + signalp_end + 1 , cursor + len(chain_cifmol.residues))))
        cursor += len(chain_cifmol.residues)
    filtered_cifmol = cifmol.residues[valid_residue_indices].extract()
    return cast("dict", filtered_cifmol.to_dict()) if len(filtered_cifmol.residues) > 0 else None

def filter_a3m(
    max_msa_depth: int = 16_384,
) -> Callable:
    """Filter instruction to select entries by resolution and date."""
    def worker(
        residue_container: FeatureContainer,
        chain_container: FeatureContainer,
    ) -> tuple[FeatureContainer, FeatureContainer]:
        msa_depth = residue_container["sequences"].shape[1]

        # remove database, database_id, rep_id
        chain_container_dict = chain_container.to_dict()
        species = chain_container_dict["nodes"]["species"]
        chain_container_dict = FeatureContainer.from_dict(FeatureContainerDict({"nodes" : {"species": species}, "edges": {}}))

        if msa_depth < max_msa_depth:
            return (
                residue_container,
                chain_container,
            )
        sequences = residue_container["sequences"]
        deletions = residue_container["deletions"]
        gap_fraction = np.mean(
            residue_container["sequences"].value == "31", # gap character in a3m
            axis=1,
        )
        sorted_indices = np.argsort(gap_fraction)
        selected_indices = sorted_indices[:max_msa_depth]
        sequences = sequences[:, selected_indices]
        sequences = NodeFeature(sequences.value.astype(np.uint8))
        deletions = deletions[:, selected_indices]
        species = species["value"][selected_indices]
        species = NodeFeature(np.array(species))
        residue_container = FeatureContainer(
            features ={
                "query_sequence": residue_container["query_sequence"],
                "sequences": sequences,
                "deletions": deletions,
                "deletion_mean": residue_container["deletion_mean"],
                "profile": residue_container["profile"],
            },
        )
        chain_container = FeatureContainer(
            features = {
                "species": species,
            },
        )

        return (
            residue_container,
            chain_container,
        )
    return worker

