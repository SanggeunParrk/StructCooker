import fnmatch
import os
import re
import subprocess
import time
from collections.abc import Callable
from pathlib import Path
from typing import cast

import kalign
import numpy as np
from biomol.core.index import IndexTable

from pipelines.cifmol import CIFMol
from pipelines.utils.io import load_bytes, load_raw_data


def _is_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def _run_command(
    command: list[str],
    *,
    env: dict[str, str] | None = None,
) -> None:
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    if result.returncode != 0:
        cmd = " ".join(command)
        msg = (
            f"Command failed ({result.returncode}): {cmd}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        raise RuntimeError(msg)


def load_a3m_list(
    data_dir: Path,
    output_dir: Path,
    pattern: str = "P*.a3m",
    output_pattern: str = ".hhm",
) -> list[dict[str, Path]]:
    """Scan a directory recursively for A3M files matching a pattern and return a list of input/output path pairs."""
    result = []

    def _scan(dir_path: Path) -> None:
        with os.scandir(dir_path) as it:
            for entry in it:
                if entry.is_dir(follow_symlinks=False):
                    _scan(Path(entry.path))
                elif fnmatch.fnmatch(entry.name, pattern):
                    result.append(
                        {
                            "input_a3m_path": Path(entry.path),
                            "output_path": output_dir
                            / f"{Path(entry.name).stem}{output_pattern}",
                        },
                    )

    _scan(data_dir)
    return result


def run_hhmake(input_a3m_path: Path, output_path: Path | None) -> str:
    """Run hhmake to convert an A3M file to an HMM file."""
    if output_path is None:
        output_path = input_a3m_path.with_suffix(".hhm")
    command = [
        "hhmake",
        "-i",
        str(input_a3m_path),
        "-o",
        str(output_path),
    ]
    try:
        _run_command(command)
        return "hhm file created at: " + str(output_path)
    except Exception as e:
        msg = f"Error running hhmake for {input_a3m_path}: {e}"
        raise RuntimeError(msg) from e


def _hhsearch_command(
    db_template: Path,
    cpu: int,
    mem: int,
    in_msa: Path,
    out_hhr: Path,
    out_atab: Path,
) -> list[str]:
    return [
        "hhsearch",
        "-b",
        "50",
        "-B",
        "500",
        "-z",
        "50",
        "-Z",
        "500",
        "-mact",
        "0.05",
        "-cpu",
        str(cpu),
        "-maxmem",
        str(mem),
        "-aliw",
        "100000",
        "-e",
        "100",
        "-p",
        "5.0",
        "-d",
        str(db_template),
        "-i",
        str(in_msa),
        "-o",
        str(out_hhr),
        "-atab",
        str(out_atab),
        "-v",
        "0",
    ]


def run_hhsearch(
    msa_path: Path,
    hhr_path: Path,
    *,
    cpu: int = 4,
    mem: int = 20,
    db_template: Path,
    hhsuite_bin_dir: Path,
) -> str:
    """Run HHsearch for one MSA directory."""
    db_template = Path(db_template)
    hhsuite_bin_dir = Path(hhsuite_bin_dir)

    hhsuite_env = os.environ.copy()
    hhsuite_env["HHLIB"] = str(hhsuite_bin_dir)
    hhsuite_env["PATH"] = f"{hhsuite_bin_dir}:{hhsuite_env.get('PATH', '')}"

    if not _is_nonempty(msa_path):
        msg = f"MSA file does not exist or is empty: {msa_path}"
        raise FileNotFoundError(msg)

    if _is_nonempty(hhr_path):
        return f"Skip {hhr_path.name} (already exists and is non-empty)"

    _run_command(
        _hhsearch_command(
            db_template=db_template,
            cpu=cpu,
            mem=mem,
            in_msa=msa_path,
            out_hhr=hhr_path,
            out_atab=hhr_path.with_suffix(".atab"),
        ),
        env=hhsuite_env,
    )
    print(f"HHsearch completed for {msa_path}, output saved to {hhr_path}")
    return f"Done {hhr_path.name}"


def run_hmmbuild(input_a3m_path: Path, hmm_path: Path | None) -> str:
    # run hmmbuild to convert a3m to hmm
    if hmm_path is None:
        hmm_path = input_a3m_path.with_suffix(".hmm")
    command = [
        "hmmbuild",
        str(hmm_path),
        str(input_a3m_path),
    ]
    if hmm_path.exists() and hmm_path.stat().st_size > 0:
        print(
            f"HMM file {hmm_path} already exists and is non-empty. Skipping hmmbuild for {input_a3m_path}.",
        )
        return f"Skip {hmm_path.name} (already exists and is non-empty)"
    try:
        _run_command(command)
        return "hmm file created at: " + str(hmm_path)
    except Exception as e:
        msg = f"Error running hmmbuild for {input_a3m_path}: {e}"
        raise RuntimeError(msg) from e


def run_hmmsearch(
    output_dir: Path,
    hmm_path: Path,
    fasta_path: Path,
) -> str:
    if not _is_nonempty(hmm_path):
        msg = f"HMM file does not exist or is empty: {hmm_path}"
        raise FileNotFoundError(msg)
    if not _is_nonempty(fasta_path):
        msg = f"FASTA file does not exist or is empty: {fasta_path}"
        raise FileNotFoundError(msg)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{hmm_path.stem}.out"
    if output_path.exists() and output_path.stat().st_size > 0:
        print(
            f"Output file {output_path} already exists and is non-empty. Skipping hmmsearch for {hmm_path}.",
        )
        return f"Skip {output_path.name} (already exists and is non-empty)"
    command = [
        "hmmsearch",
        "--noali",
        "--F1",
        "0.1",
        "--F2",
        "0.1",
        "--F3",
        "0.1",
        "-E",
        "100",
        "--incE",
        "100",
        "--domE",
        "100",
        "--incdomE",
        "100",
        "-o",
        str(output_path),
        str(hmm_path),
        str(fasta_path),
    ]
    try:
        _run_command(command)
        return f"hmmsearch completed for {hmm_path} against {fasta_path}"
    except Exception as e:
        msg = f"Error running hmmsearch for {hmm_path}: {e}"
        raise RuntimeError(msg) from e


def remove_lower_from_a3m(input_a3m_path: Path, output_path: Path | None) -> str:
    if output_path is None:
        output_path = input_a3m_path.with_suffix(".no_lower.a3m")
    try:
        with input_a3m_path.open("r") as infile, output_path.open("w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    outfile.write(line)
                else:
                    # remove all lowercase letters from the sequence lines
                    line = "".join(c for c in line if not c.islower())
                    outfile.write(line)
        return (
            f"Lowercase letters removed from {input_a3m_path}, saved to {output_path}"
        )
    except Exception as e:
        msg = f"Error processing {input_a3m_path} to remove lowercase letters: {e}"
        raise RuntimeError(msg) from e


def parse_hmm_query_mapping(hmm_path: Path) -> tuple[dict[int, int], dict[int, int]]:
    """
    Parse an HMMER3 .hmm file with 'MAP yes' and build:
      1) query_to_hmm: query/MSA column index -> HMM index
      2) hmm_to_query: HMM index -> query/MSA column index

    Returns
    -------
    query_to_hmm, hmm_to_query
    """
    query_to_hmm: dict[int, int] = {}
    hmm_to_query: dict[int, int] = {}

    with hmm_path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.rstrip()

            # node line pattern:
            # <hmm_idx> <20 emission numbers> <map_idx> <consensus> <RF> <MM> <CS>
            # example:
            # 1  2.54915 ... 3.68851      3 l - - -
            #
            # capture:
            #   group(1) = hmm_idx
            #   group(2) = map_idx (query/MSA column)
            #   group(3) = consensus residue
            m = re.match(
                r"^\s*(\d+)"  # HMM index
                r"(?:\s+\S+){20}"  # 20 match emission scores
                r"\s+(\d+)\s+([A-Za-z])"  # MAP index, consensus residue
                r"(?:\s+\S+){3}\s*$",  # RF MM CS
                s,
            )
            if m:
                hmm_idx = int(m.group(1)) - 1  # convert to 0-based index
                query_idx = int(m.group(2)) - 1

                hmm_to_query[hmm_idx] = query_idx
                query_to_hmm[query_idx] = hmm_idx

    return query_to_hmm, hmm_to_query


def extract_sequences(
    hmmsearch_output_path: Path,
    seqid2earliest_date: dict[str, time.struct_time],
    seqid2seq: dict[str, str],
) -> tuple[list[str], time.struct_time, str]:
    """Extract template chain IDs from HMMER3 hmmsearch output."""
    query_seq_id = hmmsearch_output_path.stem
    query_seq = seqid2seq.get(query_seq_id)
    if query_seq is None:
        msg = f"Query sequence not found for query sequence ID '{query_seq_id}'."
        raise KeyError(msg)
    earliest_query_date = seqid2earliest_date.get(query_seq_id)
    if earliest_query_date is None:
        msg = f"Earliest query date not found for query sequence ID '{query_seq_id}'."
        raise KeyError(msg)
    with hmmsearch_output_path.open("r", encoding="utf-8") as f:
        hmmsearch_output = f.read()
    pattern = re.compile(
        r"""
        ^\s*
        [0-9.eE+-]+      # full seq E-value
        \s+[0-9.]+       # score
        \s+[0-9.]+       # bias
        \s+[0-9.eE+-]+   # best domain E-value
        \s+[0-9.]+       # score
        \s+[0-9.]+       # bias
        \s+[0-9.]+       # exp
        \s+\d+           # N
        \s+(\S+)         # <-- Sequence (capture)
        \s+\|            # description 시작 (| 로 보장)
        """,
        re.MULTILINE | re.VERBOSE,
    )

    ids = pattern.findall(hmmsearch_output)

    chain_ids = []
    for _id in ids:
        # id format: <pdb_id>_<chain_id>
        parts = _id.split("_")
        chain_id = "_".join(
            parts[:2],
        )  # keep only the first two parts to get <pdb_id>_<chain_id>
        if chain_id not in chain_ids:
            chain_ids.append(chain_id)

    return chain_ids, earliest_query_date, query_seq


def run_kalign(
    query_seq: str,
    template_seq: str,
) -> tuple[str, str, float]:
    """Run Kalign and parse the alignment to extract aligned sequences."""
    sequences = [query_seq, template_seq]

    # Default mode — consistency anchors + VSM (best general-purpose)
    aligned = kalign.align(sequences)
    aligned_query_seq, aligned_template_seq = aligned[0], aligned[1]
    cover = 0
    for q, t in zip(aligned_query_seq, aligned_template_seq, strict=True):
        if q != "-" and t != "-":
            cover += 1
    coverage = cover / len(query_seq)
    return aligned_query_seq, aligned_template_seq, coverage


def filter_and_align_template_chain_ids(
    date_cutoff: str | None = None,
    day_diff_cutoff: int = 60,
    min_seq_len: int = 10,
    min_query_coverage: float = 0.1,
    max_query_coverage: float = 0.95,
    topk: int = 20,
) -> Callable[..., dict[str, tuple[str, str]]]:
    """Return a function that filters and aligns template chain IDs based on metadata and query date."""
    date_cutoff_time = (
        time.strptime(date_cutoff, "%Y-%m-%d") if date_cutoff is not None else None
    )

    def _worker(
        query_seq: str,
        earliest_query_date: time.struct_time,
        template_chain_ids: list[str],
        metadata_dict: dict[str, dict],  # already filtered out signal peptide
    ) -> dict[str, tuple[str, str]]:
        """Worker function."""
        align_results = {}
        for chain_id in template_chain_ids:
            metadata = metadata_dict.get(chain_id)
            if metadata is None:
                msg = f"Metadata not found for chain ID {chain_id}."
                raise KeyError(msg)
            deposit_date = metadata.get("deposition_date")
            seq = metadata.get("sequence")
            if deposit_date is None or seq is None:
                msg = f"Deposit date or sequence not found in metadata for chain ID {chain_id}."
                raise KeyError(msg)

            # 1. Filter by date cutoff
            if date_cutoff_time is not None and deposit_date > date_cutoff_time:
                continue
            # 2. Filter by day difference cutoff
            day_diff = (
                time.mktime(deposit_date) - time.mktime(earliest_query_date)
            ) / (24 * 3600)
            if day_diff < day_diff_cutoff:
                continue

            # 3. Filter by sequence length
            if len(seq) < min_seq_len:
                continue

            # 4. Filter by query coverage
            aligned_query_seq, aligned_template_seq, coverage = run_kalign(
                query_seq=query_seq,
                template_seq=seq,
            )
            if not (min_query_coverage <= coverage <= max_query_coverage):
                continue

            align_results[chain_id] = (aligned_query_seq, aligned_template_seq)
            if len(align_results) >= topk:
                break
        return align_results

    return _worker


def load_cifmol(db_path: Path, pdb_id: str, chain_id: str) -> CIFMol:
    """Load the most valuable CIFMol from LMDB by cif_id."""
    value = load_raw_data(pdb_id, db_path)

    if value is None:
        msg = f"Key '{pdb_id}' not found in LMDB database at '{db_path}'."
        raise KeyError(msg)

    value = load_bytes(value)
    max_occup_sum = -999
    best_cifmol = None

    value, metadata = value["assembly_dict"], value["metadata_dict"]

    for cif_key, _item in value.items():
        assembly_id, model_id, alt_id = cif_key.split("_")

        md = dict(metadata)
        md["assembly_id"] = assembly_id
        md["model_id"] = model_id
        md["alt_id"] = alt_id

        item = dict(_item)
        item["metadata"] = md
        item = cast("BioMolDict", item)

        cifmol = CIFMol.from_dict(item)

        chain_ids = cifmol.chains.chain_id.value
        chain_ids = {
            chain_id.split("_")[0] for chain_id in chain_ids
        }  # ignore _1, _2 etc.
        if chain_id not in chain_ids:
            continue

        occup = cifmol.atoms.occupancy.value
        # nan -> 0
        occup = np.nan_to_num(occup, nan=0.0)
        occup_sum = sum(occup) if occup is not None else 0
        if occup_sum > max_occup_sum:
            max_occup_sum = occup_sum
            best_cifmol = cifmol

    if best_cifmol is None:
        msg = (
            f"No valid CIFMol found for key '{pdb_id}' in LMDB database at '{db_path}'."
        )
        raise KeyError(msg)

    return best_cifmol


def extract_backbone_indices_from_cifmol(
    cifmol: CIFMol,
) -> np.ndarray:
    """Extract backbone atom indices (N, CA, C, CB) for each residue in the CIFMol. If CB is missing, use CA coordinates for CB."""
    residue_num = len(cifmol.residues)
    backbone_atom_indices = np.full((residue_num, 4), -1, dtype=int)
    full_xyz = cifmol.atoms.xyz

    def find_atom_index(target_xyz: np.ndarray) -> int:
        if target_xyz.size == 0:
            return -1
        matched = np.where(np.all(full_xyz == target_xyz, axis=1))[0]
        return matched[0] if matched.size > 0 else -1

    # To handle various cases of missing atoms, I think using for loop is necessary here instead of vectorized operations. We can optimize later if needed.
    for ii, residue in enumerate(cifmol.residues):
        atom_ids = residue.atoms.id.value
        xyz = residue.atoms.xyz.value

        coords = {
            atom_name: xyz[atom_ids == atom_name]
            for atom_name in ("N", "CA", "C", "CB")
        }

        backbone_atom_indices[ii, 0] = find_atom_index(coords["N"])
        backbone_atom_indices[ii, 1] = find_atom_index(coords["CA"])
        backbone_atom_indices[ii, 2] = find_atom_index(coords["C"])
        backbone_atom_indices[ii, 3] = find_atom_index(coords["CB"])

        if backbone_atom_indices[ii, 3] == -1:
            backbone_atom_indices[ii, 3] = backbone_atom_indices[ii, 1]

    return backbone_atom_indices


def to_template_mol(
    cifmol: CIFMol,
    align_result: tuple[str, str],
) -> dict:
    """Convert a CIFMol to a template mol by applying the alignment result."""
    backbone_indices = extract_backbone_indices_from_cifmol(cifmol)
    query, target = align_result

    q = np.frombuffer(query.encode(), dtype="S1")
    t = np.frombuffer(target.encode(), dtype="S1")

    q_mask = q != b"-"
    t_mask = t != b"-"

    q_idx = np.cumsum(q_mask) - 1
    t_idx = np.cumsum(t_mask) - 1

    valid = q_mask & t_mask
    q_idx = q_idx[valid]
    t_idx = t_idx[valid]

    query_seq_len = int(q_mask.sum())
    template_indices = np.full(
        (query_seq_len, 4),
        -1,
        dtype=backbone_indices.dtype,
    )
    template_indices[q_idx] = backbone_indices[t_idx]
    flattened_template_indices = template_indices.flatten()
    valid = flattened_template_indices != -1

    def _take_atom(arr: np.ndarray) -> np.ndarray:
        shape = (query_seq_len * 4, *arr.shape[1:])

        if np.issubdtype(arr.dtype, np.floating):
            fill_value = np.nan
            output = np.full(shape, fill_value, dtype=arr.dtype)
        elif np.issubdtype(arr.dtype, np.integer):
            output = np.zeros(shape, dtype=arr.dtype)
        elif np.issubdtype(arr.dtype, np.str_) or np.issubdtype(arr.dtype, np.bytes_):
            output = np.full(shape, "", dtype=arr.dtype)
        else:
            output = np.full(shape, None, dtype=object)
            arr = arr.astype(object)

        output[valid] = np.take(arr, flattened_template_indices[valid], axis=0)
        return output

    def _take_residue(arr: np.ndarray) -> np.ndarray:
        shape = (query_seq_len, *arr.shape[1:])

        if np.issubdtype(arr.dtype, np.floating):
            fill_value = np.nan
            output = np.full(shape, fill_value, dtype=arr.dtype)
        elif np.issubdtype(arr.dtype, np.integer):
            output = np.zeros(shape, dtype=arr.dtype)
        elif np.issubdtype(arr.dtype, np.str_) or np.issubdtype(arr.dtype, np.bytes_):
            output = np.full(shape, "", dtype=arr.dtype)
        else:
            output = np.full(shape, None, dtype=object)
            arr = arr.astype(object)

        output[q_idx] = np.take(arr, t_idx, axis=0)
        return output

    atom_id = _take_atom(cifmol.atoms.id.value)
    atom_xyz = _take_atom(cifmol.atoms.xyz.value)
    b_factor = _take_atom(cifmol.atoms.b_factor.value)
    occupancy = _take_atom(cifmol.atoms.occupancy.value)

    atom_dict = {
        "nodes": {
            "id": {"value": atom_id},
            "xyz": {"value": atom_xyz},
            "b_factor": {"value": b_factor},
            "occupancy": {"value": occupancy},
        },
        "edges": {},
    }

    one_letter_code_can = _take_residue(
        cifmol.residues.one_letter_code_can.value,
    )
    one_letter_code = _take_residue(cifmol.residues.one_letter_code.value)
    cif_idx = _take_residue(cifmol.residues.cif_idx.value)
    auth_idx = _take_residue(cifmol.residues.auth_idx.value)
    chem_comp_id = _take_residue(cifmol.residues.chem_comp_id.value)
    hetero = _take_residue(cifmol.residues.hetero.value)

    residue_dict = {
        "nodes": {
            "one_letter_code_can": {"value": one_letter_code_can},
            "one_letter_code": {"value": one_letter_code},
            "cif_idx": {"value": cif_idx},
            "auth_idx": {"value": auth_idx},
            "chem_comp_id": {"value": chem_comp_id},
            "hetero": {"value": hetero},
        },
        "edges": {},
    }

    entity_id = cifmol.chains.entity_id.value
    entity_type = cifmol.chains.entity_type.value
    chain_id = cifmol.chains.chain_id.value
    auth_asym_id = cifmol.chains.auth_asym_id.value
    chain_dict = {
        "nodes": {
            "entity_id": {"value": entity_id},
            "entity_type": {"value": entity_type},
            "chain_id": {"value": chain_id},
            "auth_asym_id": {"value": auth_asym_id},
        },
        "edges": {},
    }
    index_table = IndexTable.from_parents(
        atom_to_res=np.array(
            [res_idx for res_idx in range(query_seq_len) for _ in range(4)],
            dtype=int,
        ),
        res_to_chain=np.zeros(query_seq_len, dtype=int),
        n_chain=len(cifmol.chains),
    )
    metadata = cifmol.metadata

    new_dict = {
        "atoms": atom_dict,
        "residues": residue_dict,
        "chains": chain_dict,
        "index_table": index_table.to_dict(),
        "metadata": metadata,
    }
    return new_dict


def find_first(prefix: str, arr: np.ndarray) -> str | None:
    """Find the first string in arr that starts with the given prefix."""
    mask = np.char.startswith(arr, prefix)
    return arr[mask][0] if np.any(mask) else None


def load_templates(
    cif_db_path: Path,
    align_results: dict[str, tuple[str, str]],
) -> dict:
    """Load CIFMol from LMDB by cif_id."""
    if len(align_results) == 0:
        return {}
    template_mols = {}
    for full_id, align_result in align_results.items():
        pdb_id, chain_id = full_id.split("_")
        cifmol = load_cifmol(cif_db_path, pdb_id.lower(), chain_id)
        chain_id = find_first(f"{chain_id}_", cifmol.chains.chain_id.value)
        try:
            cifmol = cifmol.chains[cifmol.chains.chain_id == chain_id].extract()
        except Exception as e:
            breakpoint()
            print(f"Error occurred while extracting chain for {full_id}: {e}")
            continue
        template_mol = to_template_mol(cifmol, align_result)
        template_mols[full_id] = template_mol
    return template_mols
