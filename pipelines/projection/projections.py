from pathlib import Path


def write_fasta(data: dict[str, dict[str, dict[str, str]]], output_path: Path) -> None:
    """Write fasta dictionary to a fasta file."""
    merged_fastas = []
    for fasta_dict in data.values():
        merged_fastas.extend(fasta_dict["fasta"].values())
    merged_fasta_path = output_path / "merged.fasta"
    with merged_fasta_path.open("w") as f:
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
    data: dict[str, set[tuple[str, str]]],
    output_path: Path,
) -> None:
    """Write filtered sequence IDs to a tab-delimited file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        for src_seq_id, dst_seq_id in data["filtered_seq_ids"]:
            f.write(f"{src_seq_id}\t{dst_seq_id}\n")
