from collections.abc import Callable
from datetime import date
from typing import cast

import networkx as nx
import numpy as np

from structcooker.instructions.transforms.sequence import extract_sequence_from_cifmol
from structcooker.mols import CIFMol, CIFMolAttached
from structcooker.utils.mapping import mol_type_map, polymer_cluster_types


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
    gap_character: int = 31,
) -> Callable:
    """Filter instruction to select entries by resolution and date."""

    def worker(
        sequences: dict[str, np.ndarray],
        headers: dict[str, np.ndarray],
    ) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        aligned_sequences, deletions = (
            sequences["aligned_sequences"],
            sequences["deletions"],
        )
        msa_depth = aligned_sequences.shape[0]

        if msa_depth < max_msa_depth:
            return (
                sequences,
                headers,
            )

        gap_fraction = np.mean(
            aligned_sequences == gap_character,  # gap character in a3m
            axis=1,
        )
        sorted_indices = np.argsort(gap_fraction)
        selected_indices = sorted_indices[:max_msa_depth]
        aligned_sequences = aligned_sequences[selected_indices, :]
        deletions = deletions[selected_indices, :]

        sequences = {
            "query_sequence": sequences["query_sequence"],
            "aligned_sequences": aligned_sequences,
            "deletions": deletions,
            "deletion_mean": sequences["deletion_mean"],
            "profile": sequences["profile"],
        }

        headers = {
            "database": headers["database"][selected_indices],
            "database_id": headers["database_id"][selected_indices],
            "species": headers["species"][selected_indices],
            "rep_id": headers["rep_id"][selected_indices],
        }

        return sequences, headers

    return worker


def filter_cifmol_by_token_count(
    cifmol: CIFMol | None,
    min_token_count: int = 1,
    max_token_count: int = 512,
) -> CIFMol | None:
    """Filter instruction to remove entries with sequence token length above cutoff."""
    if cifmol is None:
        return None
    token_count = len(cifmol.residues)
    if token_count < min_token_count or token_count > max_token_count:
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
    polymer_chain_count = np.isin(
        cifmol.chains.entity_type.value,
        list(polymer_entity_types),
    ).sum()
    if polymer_chain_count > max_polymer_chain_count:
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


def filter_min_unresolved_residues(
    cifmol: CIFMol | None,
    min_unresolved_residues: int = 40,
) -> CIFMol | None:
    """Filter instruction to remove entries with unresolved residues."""
    if cifmol is None:
        return None
    valid_xyz = cifmol.atoms.xyz.value
    valid_mask = np.isfinite(valid_xyz).all(axis=1)
    valid_residue_mask = np.zeros(len(cifmol.residues), dtype=bool)
    valid_residue_mask[cifmol.index_table.atom_to_res[valid_mask]] = True
    unresolved_residue_count = np.sum(~valid_residue_mask)
    if unresolved_residue_count < min_unresolved_residues:
        return None
    return cifmol
