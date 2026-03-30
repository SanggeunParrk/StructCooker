import numpy as np

from pipelines.cifmol import CIFMolAttached
from pipelines.constants import cluster_maps


def convert_to_type(
    cluster_rep: str,
    *,
    merge_ligands: bool = True,
) -> str:
    """Convert a molecule type to its cluster representation.

    Args:
    cluster_rep (str): The cluster representation of the molecule type. (e.g. "P")

    Returns
    -------
    str: The molecule type corresponding to the cluster representation. (e.g. "polypeptide(L)")
    """
    mol_types = cluster_maps.get(cluster_rep, "Unknown")
    if merge_ligands and mol_types in {"Branched", "Non-polymer", "Unknown"}:
        return "Ligand"
    return mol_types


def analyze_db_profile(
    cifmol_dict: dict[str, dict[str, CIFMolAttached]],
) -> tuple[set, set]:
    """Analyze the database profile of a CIFMolAttached object.

    Args:
    cifmol_dict (dict[str, dict[str, CIFMolAttached]]): A dictionary of CIFMolAttached objects to analyze.

    Returns
    -------
    dict: A dictionary containing the analysis results.
    """
    # Placeholder for actual analysis logic
    monomer_clusters: set[str] = set()
    interface_clusters: set[tuple[str, str]] = set()
    for _item in cifmol_dict.values():
        cifmol = _item["cifmol"]
        cluster_ids = cifmol.chains.cluster_id.value
        cluster_ids = {str(cid) for cid in cluster_ids}
        monomer_clusters.update(cluster_ids)
        src_indices, dst_indices = (
            cifmol.chains.contact.src_indices,
            cifmol.chains.contact.dst_indices,
        )
        for src_idx, dst_idx in zip(src_indices, dst_indices, strict=False):
            src_cluster = str(cifmol.chains.cluster_id.value[src_idx])
            dst_cluster = str(cifmol.chains.cluster_id.value[dst_idx])
            src_cluster, dst_cluster = sorted([src_cluster, dst_cluster])
            interface_clusters.add((src_cluster, dst_cluster))
    return monomer_clusters, interface_clusters


def extract_edge_node(
    cifmol_dict: dict[str, dict[str, CIFMolAttached]],
) -> tuple[dict, dict]:
    """Extract edge and node information from a CIFMolAttached object.

    Args:
    cifmol_dict (dict[str, dict[str, CIFMolAttached]]): A dictionary of CIFMolAttached objects to analyze.

    Returns
    -------
    dict: A dictionary containing the extracted edge and node information.
    """
    # Placeholder for actual extraction logic
    monomer_map: dict[str, list[str]] = {}
    interface_map: dict[tuple[str, str], list[tuple[str, str]]] = {}
    for _item in cifmol_dict.values():
        cifmol = _item["cifmol"]
        valid_chains = []
        for chain_id in cifmol.chains.chain_id.value:
            xyz = cifmol.chains[
                cifmol.chains.chain_id.value == chain_id
            ].atoms.xyz.value
            valid = np.isfinite(xyz).all(axis=1)
            if valid.any():
                valid_chains.append(chain_id)
        pdb_id, assembly_id, model_id, alt_id = (
            cifmol.id[0],
            cifmol.assembly_id,
            cifmol.model_id,
            cifmol.alt_id,
        )
        cif_id = f"{pdb_id}_{assembly_id}_{model_id}_{alt_id}"
        cluster_ids = cifmol.chains.cluster_id.value
        chains = cifmol.chains.chain_id.value
        for _cluster_id, _chain_id in zip(cluster_ids, chains, strict=False):
            if _chain_id not in valid_chains:
                continue
            cluster_id = str(_cluster_id)
            chain_id = str(_chain_id)
            if cluster_id not in monomer_map:
                monomer_map[cluster_id] = []
            monomer_map[cluster_id].append(f"{cif_id}_({chain_id})")
        src_indices, dst_indices = (
            cifmol.chains.contact.src_indices,
            cifmol.chains.contact.dst_indices,
        )
        for src_idx, dst_idx in zip(src_indices, dst_indices, strict=False):
            if (
                cifmol.chains.chain_id.value[src_idx] not in valid_chains
                or cifmol.chains.chain_id.value[dst_idx] not in valid_chains
            ):
                continue
            src_cluster = str(cifmol.chains.cluster_id.value[src_idx])
            dst_cluster = str(cifmol.chains.cluster_id.value[dst_idx])
            src_cluster, dst_cluster = sorted([src_cluster, dst_cluster])
            if (src_cluster, dst_cluster) not in interface_map:
                interface_map[(src_cluster, dst_cluster)] = []
            interface_map[(src_cluster, dst_cluster)].append(
                (
                    f"{cif_id}_({cifmol.chains.chain_id.value[src_idx]!s})",
                    f"{cif_id}_({cifmol.chains.chain_id.value[dst_idx]!s})",
                ),
            )
    return monomer_map, interface_map
