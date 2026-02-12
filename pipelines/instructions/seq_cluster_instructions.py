# ruff: noqa: N802

import os
from pathlib import Path
from typing import TypeVar

from anarci import run_anarci

InputType = TypeVar("InputType", str, int, float)
FeatureType = TypeVar("FeatureType")
NumericType = TypeVar("NumericType", int, float)


def run_mmseqs2(  # noqa: PLR0913
    fasta_path: Path,
    tmp_dir: Path,
    mmseqs2_seq_id: float = 0.3,
    mmseqs2_cov: float = 0.8,
    mmseqs2_covmode: str = "0",
    mmseqs2_clustermode: str = "1",
) -> None:
    """Run MMSeqs2 to cluster sequences."""
    print(f"Running MMSeqs2 clustering for {fasta_path}...")
    os.system(  # noqa: S605
        f"mmseqs easy-cluster {fasta_path} {tmp_dir} {tmp_dir}/tmp/ "  # noqa: S108
        f"--min-seq-id {mmseqs2_seq_id} -c {mmseqs2_cov} --cov-mode {mmseqs2_covmode} --cluster-mode {mmseqs2_clustermode}",
    )
    print(f"MMSeqs2 clustering finished for {fasta_path}")


def load_fasta(fasta_path: Path) -> dict[str, str]:
    """Load fasta file into a dictionary."""
    fasta_dict = {}
    with fasta_path.open("r") as f:
        lines = f.readlines()
    current_header = ""
    for _line in lines:
        line = _line.strip()
        if line.startswith(">"):
            current_header = line[1:]
            fasta_dict[current_header] = ""
        else:
            fasta_dict[current_header] += line
    return fasta_dict


def separate_sequences(
    tmp_dir: Path,
    seq_id_map_path: Path,
    fasta_path: Path,
    sabdab_summary_path: Path,
) -> dict[str, Path]:
    """Separate sequences into different types and write to different fasta files."""
    # 1. load seq_id_map and fasta
    seq_id_dict = {}
    with seq_id_map_path.open("r") as f:
        for _line in f:
            line = _line.strip()
            seq_id, sequence = line.split("\t")
            seq_id_dict[seq_id] = sequence

    # 2. Parse SabDab summary to get antibody sequences
    fasta_dict = load_fasta(fasta_path)
    chain_id_to_seq = {}
    auth_id_to_seq = {}  # protein only (this is for sabdab)
    for header in fasta_dict:
        sequence = fasta_dict[header]
        chain_id = header.split("|")[0].strip()
        chain_id_to_seq[chain_id] = sequence
        entity_type = header.split("|")[1].strip()
        if entity_type != "polypeptide(L)":
            continue
        auth_id = header.split("Auth:")[1].strip()
        auth_id = chain_id.split("_")[0] + "_" + auth_id
        auth_id_to_seq[auth_id] = sequence

    ab_H_seq_set, ab_L_seq_set = set(), set()
    with sabdab_summary_path.open("r") as f:
        lines = f.readlines()
    for _line in lines[1:]:
        line = _line.strip()
        items = line.split("\t")
        pdb_id, ab_H_chain, ab_L_chain = items[0], items[1], items[2]
        ab_H_chain = f"{pdb_id}_{ab_H_chain}".upper() if ab_H_chain != "NA" else None
        ab_L_chain = f"{pdb_id}_{ab_L_chain}".upper() if ab_L_chain != "NA" else None
        if ab_H_chain in auth_id_to_seq:
            ab_H_seq_set.add(auth_id_to_seq[ab_H_chain])
        if ab_L_chain in auth_id_to_seq:
            ab_L_seq_set.add(auth_id_to_seq[ab_L_chain])
    # 3. Divide sequences into protein, peptide, nucleotide etc.
    entity_types = [
        "antibody_H",
        "antibody_L",
        "protein",
        "peptide",
        "protein_D",
        "nucleotide",
        "ligand",
    ]
    entity_dict = {etype: [] for etype in entity_types}
    for seq_id, sequence in seq_id_dict.items():
        if sequence in ab_H_seq_set:
            etype = "antibody_H"
        elif sequence in ab_L_seq_set:
            etype = "antibody_L"
        else:
            entity_identifier = seq_id[0]
            match entity_identifier:
                case "P":
                    etype = "protein"
                case "Q":
                    etype = "protein_D"
                case "N":
                    etype = "nucleotide"
                case "R":
                    etype = "nucleotide"
                case "D":
                    etype = "nucleotide"
                case "L":
                    etype = "ligand"
                case "B":
                    etype = "ligand"
                case "X":
                    etype = "ligand"
                case _:
                    etype = "ligand"
            if etype in ("protein", "protein_D") and len(sequence) < 10:  # noqa: PLR2004
                etype = "peptide"
        entity_dict[etype].append(seq_id)
    fasta_tmp_dir = tmp_dir / "fasta"
    fasta_tmp_dir.mkdir(parents=True, exist_ok=True)
    fasta_path_dict = {}
    for etype, seq_id_list in entity_dict.items():
        fasta_path = fasta_tmp_dir / f"{etype}.fasta"
        fasta_path_dict[etype] = fasta_path
        with fasta_path.open("w") as f:
            for seq_id in seq_id_list:
                sequence = seq_id_dict[seq_id]
                f.write(f">{seq_id}\n{sequence}\n")
    return fasta_path_dict


chotia_map = {
    "L1": list(range(24, 34 + 1)),
    "L2": list(range(50, 56 + 1)),
    "L3": list(range(89, 97 + 1)),
    "H1": list(range(26, 32 + 1)),
    "H2": list(range(50, 65 + 1)),
    "H3": list(range(95, 102 + 1)),
}


def _extract_H3L3_sequence(cdr_type: str, sequence: str) -> str:
    """Return H3 or L3 sequence from full antibody sequence using ANARCI."""
    chotia_idx = chotia_map[cdr_type]
    result = run_anarci(sequence, scheme="chothia", ncpu=16)
    if result[1][0] is None:
        return ""
    result = result[1][0][0][0]
    output = ""
    for res in result:
        idx = res[0][0]
        seq = res[1]
        if idx in chotia_idx:
            output += seq
    # remove gaps
    output = output.replace("-", "")
    return output


def extract_H3L3_sequence(
    full_fasta: Path,
    output_fasta: Path,
    cdr_type: str,
) -> list[str]:
    """Extract H3 or L3 sequences from full antibody sequences in a fasta file."""
    # if light chain only, get L3
    Ab_chain = {}
    with full_fasta.open("r") as f:
        lines = f.readlines()

    seq_id = None
    for line in lines:
        if line.startswith(">"):
            seq_id = line.strip()[1:]
        else:
            sequence = line.strip()
            Ab_chain[seq_id] = sequence

    H3L3_chain = {}
    failed_seq_ids = []
    for seq_id, sequence in Ab_chain.items():
        wo_unknown_sequence = sequence.replace("X", "")
        H3L3_sequence = _extract_H3L3_sequence(cdr_type, wo_unknown_sequence)
        if len(H3L3_sequence) == 0:
            failed_seq_ids.append(seq_id)
            continue
        H3L3_chain[seq_id] = H3L3_sequence

    with output_fasta.open("w") as f:
        for seq_id, sequence in H3L3_chain.items():
            f.write(f">{seq_id}_{cdr_type}\n")
            f.write(f"{sequence}\n")
    return failed_seq_ids


def run_CDHIT(
    input_fasta: Path,
    output_path: Path,
    seq_id: float = 0.9,
) -> None:
    """
    Run CD-HIT to cluster sequences.

    Other parameters are set to recommended values for antibody clustering.
    """
    os.system(  # noqa: S605
        f"cd-hit -i {input_fasta} -o {output_path} -c {seq_id} -n 2 -M 256000 -d 0 -T 0 -l 1 -s 0.8 -aL 0.8 -aS 0.8 -g 1",
    )


def parse_cdhit_cluster(
    cdhit_cluster_dat: Path,
) -> dict[str, list[str]]:
    """Parse CD-HIT cluster .dat file into a dictionary mapping representative to members."""
    cluster_dict = {}
    current_rep = ""
    with cdhit_cluster_dat.open("r") as f:
        for _line in f:
            line = _line.strip()
            if line.startswith(">Cluster"):
                continue
            parts = line.split()
            seq_info = parts[2]
            seq_id = seq_info.split("...")[0]
            seq_id = seq_id.replace(">", "").replace("_H3", "").replace("_L3", "")
            if line.startswith("0"):
                current_rep = seq_id.replace("P", "A")  # representative ID P->A
                cluster_dict[current_rep] = [seq_id]
            else:
                cluster_dict[current_rep].append(seq_id)
    return cluster_dict


def antibody_cluster(
    tmp_dir: Path,
    fasta_path_dict: dict[str, Path],
) -> tuple[dict[str, list[str]], set[str]]:
    """Cluster antibody sequences using CD-HIT."""
    ab_heavy_fasta_path = fasta_path_dict["antibody_H"]
    ab_light_fasta_path = fasta_path_dict["antibody_L"]
    ab_H3_fasta_path = tmp_dir / "fasta" / "antibody_H3.fasta"
    ab_L3_fasta_path = tmp_dir / "fasta" / "antibody_L3.fasta"

    ab_seq_id_list = []
    with ab_heavy_fasta_path.open("r") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line.strip()[1:]
                ab_seq_id_list.append(seq_id)
    with ab_light_fasta_path.open("r") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line.strip()[1:]
                ab_seq_id_list.append(seq_id)

    failed_H3, failed_L3 = [], []
    failed_H3 = extract_H3L3_sequence(
        ab_heavy_fasta_path,
        ab_H3_fasta_path,
        cdr_type="H3",
    )
    failed_L3 = extract_H3L3_sequence(
        ab_light_fasta_path,
        ab_L3_fasta_path,
        cdr_type="L3",
    )
    ab_H_clustered_path = tmp_dir / "cdhit_ab_H_clustered.dat"
    ab_L_clustered_path = tmp_dir / "cdhit_ab_L_clustered.dat"
    run_CDHIT(ab_H3_fasta_path, ab_H_clustered_path)
    run_CDHIT(ab_L3_fasta_path, ab_L_clustered_path)
    ab_H_cluster_output = tmp_dir / "cdhit_ab_H_clustered.dat.clstr"
    ab_L_cluster_output = tmp_dir / "cdhit_ab_L_clustered.dat.clstr"
    ab_H_cluster_dict = parse_cdhit_cluster(ab_H_cluster_output)
    ab_L_cluster_dict = parse_cdhit_cluster(ab_L_cluster_output)
    ab_cluster_dict = {**ab_H_cluster_dict, **ab_L_cluster_dict}
    # psk 20251103 : 이유는 모르겠지만, 씹히는 애들이 있음 (CDHIT에서)
    output_ab = [x for v in ab_cluster_dict.values() for x in v]
    for ab_id in ab_seq_id_list:
        if ab_id not in output_ab:
            ab_cluster_dict[ab_id.replace("P", "A")] = [ab_id]
    failed_ids = set(failed_H3 + failed_L3)
    return ab_cluster_dict, failed_ids


def merge_cluster(
    fasta_path_dict: dict[str, Path],
    protein_cluster_dict: dict[str, list[str]],
    protein_d_cluster_dict: dict[str, list[str]],
    antibody_cluster_dict: dict[str, list[str]],
) -> dict[str, list[str]]:
    """Merge cluster dictionaries."""
    cluster_dict = {
        **protein_cluster_dict,
        **protein_d_cluster_dict,
        **antibody_cluster_dict,
    }

    for etype, fasta_path in fasta_path_dict.items():
        if etype in ("protein", "protein_D", "antibody_H", "antibody_L"):
            continue
        with fasta_path.open("r") as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith(">"):
                seq_id = line.strip()[1:]
                cluster_dict[seq_id] = [seq_id]
    return cluster_dict


def parse_mmseqs2_cluster(
    mmseqs2_cluster_tsv: Path,
) -> dict[str, list[str]]:
    """Parse MMSeqs2 cluster TSV file into a dictionary mapping representative to members."""
    cluster_dict = {}
    with mmseqs2_cluster_tsv.open("r") as f:
        for _line in f:
            line = _line.strip()
            representative, member = line.split("\t")
            if representative not in cluster_dict:
                cluster_dict[representative] = []
            cluster_dict[representative].append(member)
    return cluster_dict


def recover_failed_ids(
    failed_ids: set[str],
    seq_id_map_path: Path,
    fasta_path_dict: dict[str, Path],
) -> None:
    """Recover failed IDs by assigning them to their own clusters."""
    seq_id_dict = {}
    with seq_id_map_path.open("r") as f:
        for _line in f:
            line = _line.strip()
            seq_id, sequence = line.split("\t")
            seq_id_dict[seq_id] = sequence
    protein_path = fasta_path_dict["protein"]
    failed_fasta = ""
    for seq_id in failed_ids:
        sequence = seq_id_dict[seq_id]
        failed_fasta += f">{seq_id}\n{sequence}\n"
    # add failed ids to protein fasta
    with protein_path.open("a") as f:
        f.write(failed_fasta)


def protein_cluster(
    tmp_dir: Path,
    seq_id_map_path: Path,
    fasta_path_dict: dict[str, Path],
    failed_ids: set[str],
    mmseqs2_seq_id: float = 0.3,
    mmseqs2_cov: float = 0.8,
    mmseqs2_covmode: str = "0",
    mmseqs2_clustermode: str = "1",
) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Cluster protein sequences using MMSeqs2."""
    recover_failed_ids(failed_ids, seq_id_map_path, fasta_path_dict)
    protein_path = fasta_path_dict["protein"]
    protein_d_path = fasta_path_dict["protein_D"]
    mmseqs2_protein_tmp_dir = tmp_dir / "mmseqs2_protein"
    mmseqs2_protein_d_tmp_dir = tmp_dir / "mmseqs2_protein_D"
    mmseqs2_protein_tmp_dir.mkdir(parents=True, exist_ok=True)
    mmseqs2_protein_d_tmp_dir.mkdir(parents=True, exist_ok=True)
    protein_clustered_path = tmp_dir / "mmseqs2_protein_cluster.tsv"
    protein_d_clustered_path = tmp_dir / "mmseqs2_protein_D_cluster.tsv"
    run_mmseqs2(
        protein_path,
        mmseqs2_protein_tmp_dir,
        mmseqs2_seq_id,
        mmseqs2_cov,
        mmseqs2_covmode,
        mmseqs2_clustermode,
    )
    run_mmseqs2(
        protein_d_path,
        mmseqs2_protein_d_tmp_dir,
        mmseqs2_seq_id,
        mmseqs2_cov,
        mmseqs2_covmode,
        mmseqs2_clustermode,
    )
    protein_cluster_dict = parse_mmseqs2_cluster(protein_clustered_path)
    protein_d_cluster_dict = parse_mmseqs2_cluster(protein_d_clustered_path)
    return protein_cluster_dict, protein_d_cluster_dict
