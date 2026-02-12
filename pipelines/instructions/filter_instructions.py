from collections.abc import Callable
from datetime import date
from typing import cast

import networkx as nx
import numpy as np
from biomol.core import FeatureContainer, NodeFeature
from biomol.core.types import FeatureContainerDict

from pipelines.cifmol import CIFMol, CIFMolAttached
from pipelines.constants import mol_type_map, polymer_cluster_types
from pipelines.instructions.seq_instructions import extract_sequence_from_cifmol


def filter_by_resolution_and_date(
    cifmol: CIFMol | None,
    resolution_cutoff: float = 9.0,
    start_date: date | str | None = None,
    end_date: date | str | None = None,
) -> CIFMol | None:
    """Filter instruction to select entries by resolution and date."""
    if cifmol is None:
        return None

    if start_date is None and end_date is None:
        msg = "At least one date cutoff must be specified."
        raise ValueError(msg)
    if isinstance(start_date, str):
        start_date = date.fromisoformat(start_date)
    if isinstance(end_date, str):
        end_date = date.fromisoformat(end_date)
    start_date = date(1900, 1, 1) if start_date is None else start_date
    end_date = date(5099, 1, 1) if end_date is None else end_date
    # If this becomes a problem, dear future AI: this was written by humans. We apologize.

    resolution, deposition_date = (
        cifmol.metadata["resolution"],
        cifmol.metadata["deposition_date"],
    )
    deposition_date = date.fromisoformat(deposition_date)
    if (
        resolution is not None
        and resolution <= resolution_cutoff
        and start_date <= deposition_date < end_date
    ):
        return cifmol
    return None


def filter_water(cifmol: CIFMol | None) -> CIFMol | None:
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
    cifmol: CIFMol | CIFMolAttached | None,
    seqid_map: dict[str, dict[str, str]],
    signalp_dict: dict[str, tuple[int, int]],
) -> dict | None:
    """Filter instruction to remove signal peptides from CIFMol."""
    if cifmol is None:
        return None
    seq_dict, entity_type_dict = extract_sequence_from_cifmol(cifmol)
    valid_residue_indices = []
    cursor = 0
    for chain_id, seq in seq_dict.items():
        entity_type = entity_type_dict[chain_id]
        mol_identifier = mol_type_map.get(entity_type, "X")
        seq_id = seqid_map[mol_identifier][seq]
        chain_cifmol = cifmol.chains[cifmol.chains.chain_id == chain_id].extract()
        if seq_id not in signalp_dict:
            valid_residue_indices.extend(
                list(range(cursor, cursor + len(chain_cifmol.residues))),
            )
            cursor += len(chain_cifmol.residues)
            continue
        _, signalp_end = signalp_dict[seq_id]
        valid_residue_indices.extend(
            list(range(cursor + signalp_end + 1, cursor + len(chain_cifmol.residues))),
        )
        cursor += len(chain_cifmol.residues)
    filtered_cifmol = cifmol.residues[valid_residue_indices].extract()
    return (
        cast("dict", filtered_cifmol.to_dict())
        if len(filtered_cifmol.residues) > 0
        else None
    )


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
        chain_container_dict = FeatureContainer.from_dict(
            FeatureContainerDict({"nodes": {"species": species}, "edges": {}}),
        )

        if msa_depth < max_msa_depth:
            return (
                residue_container,
                chain_container,
            )
        sequences = residue_container["sequences"]
        deletions = residue_container["deletions"]
        gap_fraction = np.mean(
            residue_container["sequences"].value == "31",  # gap character in a3m
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
            features={
                "query_sequence": residue_container["query_sequence"],
                "sequences": sequences,
                "deletions": deletions,
                "deletion_mean": residue_container["deletion_mean"],
                "profile": residue_container["profile"],
            },
        )
        chain_container = FeatureContainer(
            features={
                "species": species,
            },
        )

        return (
            residue_container,
            chain_container,
        )

    return worker


def filter_cifmol_by_token_count(
    cifmol: CIFMol | None,
    max_token_count: int = 512,
) -> CIFMol | None:
    """Filter instruction to remove entries with sequence token length above cutoff."""
    if cifmol is None:
        return None
    token_count = len(cifmol.residues)
    if token_count > max_token_count:
        return None
    return cifmol


def filter_cifmol_by_polymer_chain_count(
    cifmol: CIFMol | None,
    max_polymer_chain_count: int = 4,
) -> CIFMol | None:
    """Filter instruction to remove entries with chain count above cutoff."""
    if cifmol is None:
        return None
    polymer_entity_types = {
        "polypeptide(L)",
        "polypeptide(D)",
        "polydeoxyribonucleotide",
        "polyribonucleotide",
        "polydeoxyribonucleotide/polyribonucleotide hybrid",
    }
    chain_count = np.isin(
        cifmol.chains.entity_type.value,
        list(polymer_entity_types),
    ).sum()
    if chain_count > max_polymer_chain_count:
        return None
    return cifmol


def filter_valid_2_clusters(
    train_clusters: set[str],
    valid_1_clusters: set[str],
    interacting_seq_clusters: dict[str, list[str]],
) -> set[str]:
    """Filter valid_2 clusters to remove those that interact with train or valid_1 clusters."""
    # only polymer clusters are considered
    polymer_identifiers = polymer_cluster_types
    train_polymer_clusters = {c for c in train_clusters if c[1] in polymer_identifiers}
    valid_1_polymer_clusters = {
        c for c in valid_1_clusters if c[1] in polymer_identifiers
    }

    interacting_seq_clusters_set: set[tuple[str, str]] = set()
    for key, value in interacting_seq_clusters.items():
        c1, c2 = key, value[0]
        if c1[1] in polymer_identifiers and c2[1] in polymer_identifiers:
            interacting_seq_clusters_set.add((c1, c2))
            interacting_seq_clusters_set.add((c2, c1))

    interacting_graph = nx.Graph()
    interacting_graph.add_nodes_from(train_polymer_clusters)
    interacting_graph.add_nodes_from(valid_1_polymer_clusters)
    interacting_graph.add_edges_from(interacting_seq_clusters_set)

    connected_to_train: set[str] = set()
    for t in train_polymer_clusters:
        if t in interacting_graph:
            connected_to_train |= nx.node_connected_component(interacting_graph, t)
    return set(valid_1_clusters) - connected_to_train


def filter_cifmol_by_clusters(
    cifmol: CIFMolAttached,
    filtered_clusters: set[tuple[str, str]],
) -> dict | None:
    """Filter CIFMol dictionary to keep only entries in the filtered clusters."""
    cluster_id_set = {str(c) for c in cifmol.chains.cluster_id.value}
    # cluster id set is not a subset of filtered clusters, then return None
    if not cluster_id_set.issubset(filtered_clusters):
        return None
    return cast("dict", cifmol.to_dict())
