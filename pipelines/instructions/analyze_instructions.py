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
        cluster_ids = set(str(cid) for cid in cluster_ids)
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
