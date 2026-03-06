from pathlib import Path
from typing import cast

import numpy as np
from biomol.core import NodeFeature

from pipelines.cifmol import CIFMol, CIFMolAttached
from pipelines.constants import mol_type_map


def load_tsv(
    tsv_file_path: Path,
    *,
    split_by_comma: bool = True,
) -> dict[str, list[str]]:
    """Load a TSV file and return its contents as a dictionary."""
    data: dict[str, list[str]] = {}
    with tsv_file_path.open("r", encoding="utf-8") as file:
        lines = file.readlines()
        for line in lines:
            key, value = line.strip().split("\t")
            if not split_by_comma:
                data[key] = [value]
            else:
                data[key] = value.split(",")
    return data


def reverse_dict(
    input_dict: dict[str, list[str]],
) -> dict[str, str]:
    """Reverse a dictionary mapping from str to list[str] into a dictionary mapping from str to str."""
    output_dict: dict[str, str] = {}
    for key, value_list in input_dict.items():
        for value in value_list:
            output_dict[value] = key
    return output_dict


def build_seqid_map(
    seqid2seq: dict[str, list[str]],
) -> dict[str, dict[str, str]]:
    """Build a mapping from sequence and molecule type to sequence ID."""
    seqid_map: dict[str, dict[str, str]] = {}  # first key : mol type, second key : seq
    for mol_identifier in mol_type_map.values():
        seqid_map[mol_identifier] = {}
    for seqid, seq_list in seqid2seq.items():
        seq = seq_list[0]
        mol_identifier = seqid[0]
        if mol_identifier not in seqid_map:
            msg = f"Unknown molecule identifier {mol_identifier} in seqid {seqid}."
            raise ValueError(msg)
        seqid_map[mol_identifier][seq] = seqid
    return seqid_map


def build_seq_metadata_map(
    raw_fasta_dict: dict[str, str],
    seqid2seq: dict[str, list[str]],
    seqclusters2seqids: dict[str, list[str]],
) -> dict[str, tuple[str, str]]:
    """Build a metadata map from sequence cluster ID to sequence ID."""
    seq_metadata_map: dict[str, tuple[str, str]] = {}  # cif_id -> (seq id, seq cluster)
    seqid2seqcluster: dict[str, str] = reverse_dict(seqclusters2seqids)
    seqid_map = build_seqid_map(seqid2seq)

    for header, sequence in raw_fasta_dict.items():
        mol_type = header.split("|")[1].strip()
        mol_identifier = mol_type_map.get(mol_type, "X")
        seqid = seqid_map[mol_identifier].get(sequence)
        if seqid is None:
            msg = f"Sequence ID not found for molecule type {mol_identifier} and sequence {sequence}."
            raise KeyError(msg)
        seqcluster = seqid2seqcluster.get(seqid)
        if seqcluster is None:
            msg = f"Sequence cluster not found for sequence ID {seqid}."
            raise KeyError(msg)
        cif_id = header.split("|")[0].strip()  # pdbid_chainid_altid
        seq_metadata_map[cif_id] = (seqid, seqcluster)
    return seq_metadata_map


def parse_signalp(signalp_path: Path) -> tuple[int, int] | None:
    """Parse the signalp output file and extract the sequence ids."""
    if not signalp_path.exists():
        return None
    with signalp_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    result = lines[1].split("\t")
    return int(result[3]) - 1, int(result[4]) - 1


def load_signalp(
    signalp_dir: Path | None,
) -> dict[str, list[tuple[int, int]]]:
    """Load SignalP data from a directory containing GFF3 files."""
    signalp_data = {}
    if signalp_dir is not None:
        if not signalp_dir.exists():
            msg = f"SignalP directory {signalp_dir} does not exist."
            raise FileNotFoundError(msg)
        for signalp_file in signalp_dir.glob("*.gff3"):
            seqid = signalp_file.stem
            signalp_data[seqid] = parse_signalp(signalp_file)
    return signalp_data


def attach_metadata(
    cifmol: CIFMol,
    seq_metadata_map: dict[str, tuple[str, str]],
) -> dict:
    """Attach metadata to a CIFMol object."""
    pdbid = cifmol.id[0]
    alt_id = cifmol.alt_id
    seq_id_list = []
    seq_cluster_list = []
    for full_chain_id in cifmol.chains.chain_id.value:
        chain_id = full_chain_id.split("_")[0]
        cif_key = f"{pdbid}_{chain_id}_{alt_id}"
        if cif_key not in seq_metadata_map:
            msg = f"Sequence metadata not found for CIF key {cif_key}."
            raise KeyError(msg)
        seqid, seqcluster = seq_metadata_map[cif_key]
        seq_id_list.append(seqid)
        seq_cluster_list.append(seqcluster)
    seq_id_list = np.array(seq_id_list)
    seq_cluster_list = np.array(seq_cluster_list)
    seq_id_list = NodeFeature(seq_id_list)
    seq_cluster_list = NodeFeature(seq_cluster_list)
    cifmol_dict = cifmol.to_dict()
    cifmol_dict["chains"]["nodes"]["seq_id"] = {
        "value": np.array(seq_id_list, dtype=str),
    }
    cifmol_dict["chains"]["nodes"]["cluster_id"] = {
        "value": np.array(seq_cluster_list, dtype=str),
    }
    return cast("dict", CIFMolAttached.from_dict(cifmol_dict).to_dict())


def extract_metadata(
    cifmol_dict: dict[str, dict[str, CIFMol]],
) -> dict[str, dict[str, str]]:
    """Extract resolution and date etc from a CIFMol object."""
    metadata_dict = {}
    NA_types = {
        "polydeoxyribonucleotide",
        "polyribonucleotide",
        "polydeoxyribonucleotide/polyribonucleotide hybrid",
    }
    D_types = {"polypeptide(D)"}

    for cif_key, cifmol_wrapper in cifmol_dict.items():
        cifmol = cifmol_wrapper["cifmol"]
        metadata_dict[cif_key] = {}
        resolution = cifmol.metadata.get("resolution", "NA")
        deposition_date = cifmol.metadata.get("deposition_date", "NA")
        chain_num = len(cifmol.chains)
        residue_num = len(cifmol.residues)
        atom_num = len(cifmol.atoms)
        including_NA = set(cifmol.chains.entity_type.value).intersection(NA_types)
        including_NA = "Yes" if including_NA else "No"
        including_Dform = set(cifmol.chains.entity_type.value).intersection(D_types)
        including_Dform = "Yes" if including_Dform else "No"
        metadata_dict[cif_key]["resolution"] = str(resolution)
        metadata_dict[cif_key]["deposition_date"] = str(deposition_date)
        metadata_dict[cif_key]["chain_num"] = str(chain_num)
        metadata_dict[cif_key]["residue_num"] = str(residue_num)
        metadata_dict[cif_key]["atom_num"] = str(atom_num)
        metadata_dict[cif_key]["including_NA"] = including_NA
        metadata_dict[cif_key]["including_Dform"] = including_Dform

    return metadata_dict


def classify_seq_clusters(
    raw_fasta_dict: dict[str, str],
    fasta_dict: dict[str, str],
    seqid_map: dict[str, dict[str, str]],
    seqclusters2seqids: dict[str, list[str]],
) -> set[str]:
    """Classify sequence clusters based on the provided fasta dictionary and sequence ID map."""
    seqid2seqcluster: dict[str, str] = reverse_dict(seqclusters2seqids)
    classified_clusters: set[str] = set()
    for header in fasta_dict:
        mol_type = header.split("|")[1].strip()
        mol_identifier = mol_type_map.get(mol_type, "X")
        raw_sequence = raw_fasta_dict[header]
        seqid = seqid_map[mol_identifier].get(raw_sequence)
        if seqid is None:
            msg = f"Sequence ID not found for molecule type {mol_identifier} and sequence {raw_sequence}."
            raise KeyError(msg)
        seqcluster = seqid2seqcluster.get(seqid)
        if seqcluster is None:
            msg = f"Sequence cluster not found for sequence ID {seqid}."
            raise KeyError(msg)
        classified_clusters.add(seqcluster)
    return classified_clusters

def load_fasta(
    fasta_path: Path,
) -> dict[str, str]:
    """Load a FASTA file and return its contents as a dictionary."""
    fasta_dict = {}
    with fasta_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        current_header = None
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                current_header = line[1:]  # Remove the '>' character
                fasta_dict[current_header] = ""
            elif current_header is not None:
                fasta_dict[current_header] += line
            else:
                msg = "FASTA format error: sequence data found before any header."
                raise ValueError(msg)
    return fasta_dict


def extract_protein_seqs(
    seqid2seq: dict[str, list[str]],
    remove_unknown: bool = True,
) -> list[dict]:
    """Extract protein sequences from the seqid2seq mapping."""
    protein_seqs = []
    for seqid, seqs in seqid2seq.items():
        if len(seqs) == 0:
            msg = f"No sequence found for sequence ID {seqid}."
            raise ValueError(msg)
        if len(seqs) > 1:
            msg = f"Multiple sequences found for sequence ID {seqid}."
            raise ValueError(msg)
        seq = seqs[0]
        mol_identifier = seqid[0]
        if mol_identifier == "P":
            if remove_unknown:
                is_unknown = all(aa == "X" for aa in seq)
                if is_unknown:
                    continue
            protein_seqs.append({"seqid": seqid, "sequence": seq})
    return protein_seqs
