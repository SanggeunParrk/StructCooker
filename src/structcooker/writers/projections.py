from pathlib import Path

from structcooker.constants import cluster_maps


def write_fasta(data: dict[str, dict[str, dict[str, str]]], output_path: Path) -> None:
    """Write fasta dictionary to a fasta file."""
    merged_fastas = []
    for fasta_dict in data.values():
        merged_fastas.extend(fasta_dict["fasta"].values())
    with output_path.open("w") as f:
        for fasta in merged_fastas:
            f.write(fasta)


def write_seq_id_map(data: dict[str, dict[str, str]], output_path: Path) -> None:
    """Write sequence hash map to a tab-delimited file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    seq_id_maps = data["seq_id_map"]
    # sort by seq_id (string)
    seq_id_maps = dict(
        sorted(seq_id_maps.items(), key=lambda item: item[0]),
    )
    with output_path.open("w") as f:
        for seq_id, sequence in seq_id_maps.items():
            f.write(f"{seq_id}\t{sequence}\n")


def write_seq_cluster_dict(data: dict, output_path: Path) -> None:
    """Write sequence cluster dictionary to a tab-delimited file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        for rep_seq_hash, member_list in data["cluster_dict"].items():
            members = ",".join(member_list)
            f.write(f"c{rep_seq_hash}\t{members}\n")


def write_metadata(
    data: dict[str, dict[str, dict[str, dict[str, str]]]],
    output_path: Path,
) -> None:
    """Write fasta dictionary to a fasta file."""
    metadata_keys = [
        "cif_id",
        "resolution",
        "deposition_date",
        "chain_num",
        "residue_num",
        "atom_num",
        "including_NA",
        "including_Dform",
    ]
    header = "\t".join(metadata_keys)
    lines = [header]

    for pdb_id, _metadata_dict in data.items():
        metadata_dict = _metadata_dict["metadata_dict"]
        for cif_key, value in metadata_dict.items():
            cif_id = pdb_id + "_" + cif_key
            line_elems = [cif_id]
            line_elems.extend([value.get(key, "NA") for key in metadata_keys[1:]])
            line = "\t".join(line_elems)
            lines.append(line)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        for line in lines:
            f.write(line + "\n")


def write_clusters(
    data: dict[str, dict[str, list[str]]],
    output_path: Path,
) -> None:
    """Write cluster dictionaries to pickle files."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        for rep_seq_hash, member_list in data["cluster_dict"].items():
            members = ",".join(member_list)
            f.write(f"c{rep_seq_hash}\t{members}\n")


def write_filtered_seq_ids(
    data: dict[str, dict[str, set[tuple[str, str]]]],
    output_path: Path,
) -> None:
    """Write filtered sequence IDs to a tab-delimited file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    filtered_seq_ids = set()
    for seq_id_dict in data.values():
        for src_seq_id, dst_seq_id in seq_id_dict["filtered_seq_ids"]:
            filtered_seq_ids.add((src_seq_id, dst_seq_id))
    # sort by src_seq_id, dst_seq_id
    filtered_seq_ids = sorted(filtered_seq_ids, key=lambda x: (x[0], x[1]))
    with output_path.open("w") as f:
        for src_seq_id, dst_seq_id in sorted(filtered_seq_ids):
            f.write(f"{src_seq_id}\t{dst_seq_id}\n")


def write_filtered_seq_clusters(
    data: dict[str, set[tuple[str, str]]],
    output_path: Path,
) -> None:
    """Write filtered sequence clusters to a tab-delimited file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    interacting_seq_clusters = set()
    for src_seq_id, dst_seq_id in data["interacting_seq_clusters"]:
        interacting_seq_clusters.add((src_seq_id, dst_seq_id))
    interacting_seq_clusters = sorted(
        interacting_seq_clusters,
        key=lambda x: (x[0], x[1]),
    )
    with output_path.open("w") as f:
        for src_seq_id, dst_seq_id in sorted(interacting_seq_clusters):
            f.write(f"{src_seq_id}\t{dst_seq_id}\n")


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


def write_statistics(
    data: dict[str, dict],
    output_path: Path,
) -> None:
    """Write statistics to a tab-delimited file."""
    whole_monomer_clusters = set()
    whole_interface_clusters = set()
    for inner_dict in data.values():
        monomer_clusters, interface_clusters = (
            inner_dict["monomer_clusters"],
            inner_dict["interface_clusters"],
        )
        whole_monomer_clusters.update(monomer_clusters)
        whole_interface_clusters.update(interface_clusters)
    monomer_stats = {}
    interface_stats = {}
    for cluster_rep in whole_monomer_clusters:
        mol_type = convert_to_type(cluster_rep[1])
        monomer_stats[mol_type] = monomer_stats.get(mol_type, 0) + 1
    for src_rep, dst_rep in whole_interface_clusters:
        src_type = convert_to_type(src_rep[1])
        dst_type = convert_to_type(dst_rep[1])
        type_pair = tuple(sorted([src_type, dst_type]))
        interface_stats[type_pair] = interface_stats.get(type_pair, 0) + 1
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        f.write("Monomer Cluster Statistics:\n")
        for mol_type, count in monomer_stats.items():
            f.write(f"{mol_type},{count}\n")
        f.write("\nInterface Cluster Statistics:\n")
        for type_pair, count in interface_stats.items():
            f.write(f"{type_pair[0]}-{type_pair[1]},{count}\n")


def write_each_ccd_cif(
    data: dict[str, dict[str, str]],
    output_path: Path,
) -> None:
    """Write each CCD CIF payload into its own file."""
    cif_dict = data["each_cif_lines"]
    for cif_id, cif_lines in cif_dict.items():
        cif_file_path = output_path / f"{cif_id[0]}" / f"{cif_id}.cif"
        cif_file_path.parent.mkdir(parents=True, exist_ok=True)
        with cif_file_path.open("w") as f:
            f.writelines(cif_lines)


def write_edge_node(
    data: dict[str, dict[str, dict]],
    output_path: Path,
) -> None:
    """Write edge and node information to a tab-delimited file."""
    header = "cluster1\tcluster2\tpdb_id\tassembly_id\tmodel_id\talt_id\tchain_id1\tchain_id2\n"
    monomer_map = {}
    interface_map = {}

    for value in data.values():
        _monomer_map, _interface_map = value["monomer_map"], value["interface_map"]

        for key in _monomer_map:
            if key not in monomer_map:
                monomer_map[key] = []
            monomer_map[key].extend(_monomer_map[key])
        for key in _interface_map:
            if key not in interface_map:
                interface_map[key] = []
            interface_map[key].extend(_interface_map[key])

    # remove cluster_to_chain_map entries if the cluster is in any interface cluster
    interface_clusters = {
        cluster for interface in interface_map for cluster in interface
    }
    single_clusters = set(monomer_map.keys()) - interface_clusters
    filtered_monomer_map = {key: monomer_map[key] for key in single_clusters}

    def _extract_info(chain_info: dict) -> tuple[str, str, str, str, str]:
        return (
            chain_info["pdb_id"],
            chain_info["assembly_id"],
            chain_info["model_id"],
            chain_info["alt_id"],
            chain_info["chain_id"],
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        f.write(header)
        for cluster1, chain_info_list in filtered_monomer_map.items():
            for chain_info in chain_info_list:
                pdb_id, assembly_id, model_id, alt_id, chain_id = _extract_info(
                    chain_info,
                )
                to_write = f"{cluster1}\tNone\t{pdb_id}\t{assembly_id}\t{model_id}\t{alt_id}\t{chain_id}\tNone\n"
                f.write(to_write)
        for (cluster1, cluster2), chain_id_pairs_list in interface_map.items():
            for chain_id_pairs in chain_id_pairs_list:
                pdb_id1, assembly_id1, model_id1, alt_id1, chain_id1 = _extract_info(
                    chain_id_pairs[0],
                )
                _pdb_id2, _assembly_id2, _model_id2, _alt_id2, chain_id2 = _extract_info(
                    chain_id_pairs[1],
                )
                to_write = f"{cluster1}\t{cluster2}\t{pdb_id1}\t{assembly_id1}\t{model_id1}\t{alt_id1}\t{chain_id1}\t{chain_id2}\n"
                f.write(to_write)


def unittest(
    data: dict[str, dict[str, str]],
    output_path: Path,
) -> None:
    """Write unit test results to a tab-delimited file."""
    with output_path.open("w") as f:
        for item_id, _result in data.items():
            result = _result["results"]
            if result is None:
                result = "passed"
            f.write(f"{item_id}\t{result}\n")
