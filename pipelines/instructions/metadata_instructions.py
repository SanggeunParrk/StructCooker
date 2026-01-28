from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np

from pipelines.cifmol import CIFMol, CIFMolAttached
from pipelines.instructions.seq_instructions import extract_sequence_from_cifmol

if TYPE_CHECKING:
    from biomol.core.types import BioMolDict

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
    seq2seqid: dict[str, str],
    seqid2clusterid: dict[str, str],
) -> CIFMolAttached | None:
    """Attach metadata to a CIFMol object."""
    if cifmol is None:
        return None

    seqid_list = []
    clusterid_list = []
    seq_dict = extract_sequence_from_cifmol(cifmol)
    for seq in seq_dict.values():
        seqid = seq2seqid[seq]
        seqid_list.append(seqid)
        clusterid = seqid2clusterid[seqid]
        clusterid_list.append(clusterid)
    cifmol_dict = cifmol.to_dict()
    cifmol_dict = cast("dict", cifmol_dict)
    cifmol_dict["chains"]["nodes"]["seq_id"] = {
        "value" : np.array(seqid_list, dtype=str),
    }
    cifmol_dict["chains"]["nodes"]["cluster_id"] = {
        "value" : np.array(clusterid_list, dtype=str),
    }
    cifmol_dict = cast("BioMolDict", cifmol_dict)
    return CIFMolAttached.from_dict(cifmol_dict)

