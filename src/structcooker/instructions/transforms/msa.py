import os
import re
import shutil
import string
import subprocess
from pathlib import Path
from typing import TypeVar

import numpy as np

from structcooker.utils.mapping import ResidueMapping

DEFAULT_DB_UR30 = Path(
    "/data/shared/cssb_data/db_protSeq/uniref30/2022_02/UniRef30_2022_02",
)
DEFAULT_DB_BFD = Path(
    "/data/shared/cssb_data/db_protSeq/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
)
DEFAULT_HHSUITE_BIN = Path("/software/hhsuite/build/bin")


def _is_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def _count_header_lines(path: Path) -> int:
    with path.open("r") as f:
        return sum(1 for line in f if line.startswith(">"))


def _run_command(
    command: list[str],
    *,
    env: dict[str, str] | None = None,
) -> None:
    result = subprocess.run(command, capture_output=True, text=True, env=env, check=False)
    if result.returncode != 0:
        cmd = " ".join(command)
        msg = (
            f"Command failed ({result.returncode}): {cmd}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        raise RuntimeError(msg)


def _hhblits_command(
    db_path: Path,
    cpu: int,
    mem: int,
    in_a3m: Path,
    out_a3m: Path,
    evalue: str,
) -> list[str]:
    return [
        "hhblits",
        "-o",
        "/dev/null",
        "-mact",
        "0.35",
        "-maxfilt",
        "100000000",
        "-neffmax",
        "20",
        "-cov",
        "25",
        "-cpu",
        str(cpu),
        "-nodiff",
        "-realign_max",
        "100000000",
        "-maxseq",
        "1000000",
        "-maxmem",
        str(mem),
        "-n",
        "4",
        "-d",
        str(db_path),
        "-i",
        str(in_a3m),
        "-oa3m",
        str(out_a3m),
        "-e",
        evalue,
        "-v",
        "0",
    ]


def _hhfilter_command(in_a3m: Path, out_a3m: Path, coverage: int) -> list[str]:
    return [
        "hhfilter",
        "-maxseq",
        "100000",
        "-id",
        "90",
        "-cov",
        str(coverage),
        "-i",
        str(in_a3m),
        "-o",
        str(out_a3m),
    ]

def make_input_fasta(
    seqid: str,
    sequence: str,
    output_dir: Path,
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    out_dir = output_dir / seqid[:4] / seqid
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = out_dir / f"{seqid}.fasta"
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with fasta_path.open("w") as f:
        f.write(f">{seqid}\n{sequence}\n")
    return fasta_path, out_dir

def run_signalp(  # noqa: PLR0913
    input_fasta: Path,
    out_dir: Path,
    *,
    signalp_mode: str = "fast",
) -> Path:
    """Run SignalP + HHblits/HHfilter for one FASTA."""
    out_dir.mkdir(parents=True, exist_ok=True)
    signalp_dir = out_dir / "signalp"
    signalp_dir.mkdir(parents=True, exist_ok=True)

    _run_command(
        [
            "signalp6",
            "--fastafile",
            str(input_fasta),
            "--organism",
            "other",
            "--output_dir",
            str(signalp_dir),
            "--format",
            "none",
            "--mode",
            signalp_mode,
        ],
    )

    trim_fasta = signalp_dir / "processed_entries.fasta"
    return trim_fasta if _is_nonempty(trim_fasta) else input_fasta

def run_msa_search(  # noqa: PLR0913
    input_fasta: Path,
    out_dir: Path,
    *,
    cpu: int = 4,
    mem: int = 20,
    db_ur30: Path,
    db_bfd: Path,
    hhsuite_bin_dir: Path,
) -> str:

    hhsuite_env = os.environ.copy()
    hhsuite_env["HHLIB"] = str(hhsuite_bin_dir)
    hhsuite_env["PATH"] = f"{hhsuite_bin_dir}:{hhsuite_env.get('PATH', '')}"

    hhblits_dir = out_dir / "hhblits"
    hhblits_dir.mkdir(parents=True, exist_ok=True)
    msa0_file = out_dir / "t000_msa0.a3m"

    if not _is_nonempty(msa0_file):
        prev_a3m = input_fasta
        for evalue in ("1e-10", "1e-6", "1e-3"):
            a3m_file = hhblits_dir / f"t000_.{evalue}.a3m"
            if not _is_nonempty(a3m_file):
                _run_command(
                    _hhblits_command(
                        db_path=db_ur30,
                        cpu=cpu,
                        mem=mem,
                        in_a3m=prev_a3m,
                        out_a3m=a3m_file,
                        evalue=evalue,
                    ),
                    env=hhsuite_env,
                )

            id90cov75_file = hhblits_dir / f"t000_.{evalue}.id90cov75.a3m"
            id90cov50_file = hhblits_dir / f"t000_.{evalue}.id90cov50.a3m"
            _run_command(
                _hhfilter_command(
                    in_a3m=a3m_file,
                    out_a3m=id90cov75_file,
                    coverage=75,
                ),
                env=hhsuite_env,
            )
            _run_command(
                _hhfilter_command(
                    in_a3m=a3m_file,
                    out_a3m=id90cov50_file,
                    coverage=50,
                ),
                env=hhsuite_env,
            )
            prev_a3m = id90cov50_file

            n75 = _count_header_lines(id90cov75_file)
            n50 = _count_header_lines(id90cov50_file)

            if n75 > 2000:  # noqa: PLR2004
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(id90cov75_file, msa0_file)
                    break
            elif n50 > 4000:  # noqa: PLR2004
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(id90cov50_file, msa0_file)
                    break

        if not _is_nonempty(msa0_file):
            evalue = "1e-3"
            bfd_a3m_file = hhblits_dir / f"t000_.{evalue}.bfd.a3m"
            if not _is_nonempty(bfd_a3m_file):
                _run_command(
                    _hhblits_command(
                        db_path=db_bfd,
                        cpu=cpu,
                        mem=mem,
                        in_a3m=prev_a3m,
                        out_a3m=bfd_a3m_file,
                        evalue=evalue,
                    ),
                    env=hhsuite_env,
                )

            bfd_id90cov75_file = hhblits_dir / f"t000_.{evalue}.bfd.id90cov75.a3m"
            bfd_id90cov50_file = hhblits_dir / f"t000_.{evalue}.bfd.id90cov50.a3m"
            _run_command(
                _hhfilter_command(
                    in_a3m=bfd_a3m_file,
                    out_a3m=bfd_id90cov75_file,
                    coverage=75,
                ),
                env=hhsuite_env,
            )
            _run_command(
                _hhfilter_command(
                    in_a3m=bfd_a3m_file,
                    out_a3m=bfd_id90cov50_file,
                    coverage=50,
                ),
                env=hhsuite_env,
            )
            prev_a3m = bfd_id90cov50_file

            n75 = _count_header_lines(bfd_id90cov75_file)
            n50 = _count_header_lines(bfd_id90cov50_file)
            if n75 > 2000:  # noqa: PLR2004
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(bfd_id90cov75_file, msa0_file)
            elif n50 > 4000:  # noqa: PLR2004
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(bfd_id90cov50_file, msa0_file)

        if not _is_nonempty(msa0_file):
            shutil.copyfile(prev_a3m, msa0_file)

    return f"Done {input_fasta}"



# ---- merged from a3m.py ----

InputType = TypeVar("InputType", str, int, float)
FeatureType = TypeVar("FeatureType")
NumericType = TypeVar("NumericType", int, float)


def parse_sequence(
    raw_sequences: list[str],
    a3m_type: str | None = "protein",
) -> dict[str, np.ndarray]:
    """Parse a sequence string into a list of residue symbols."""
    table = str.maketrans(dict.fromkeys(string.ascii_lowercase))

    residue_mapping = ResidueMapping()
    max_idx = residue_mapping.MAX_INDEX

    if a3m_type is None:
        a3m_type = "protein"

    match a3m_type.lower():
        case "protein":
            mapping_view = residue_mapping.protein
        case "rna":
            mapping_view = residue_mapping.rna
        case _:
            msg = f"Unsupported a3m_type: {a3m_type}"
            raise ValueError(msg)

    query_sequence = raw_sequences[0]
    length = len(query_sequence)

    sequences = []
    deletions = []

    for raw_sequence in raw_sequences:
        lower_case = np.array(
            [0 if c.isupper() or c == "-" else 1 for c in raw_sequence],
        )
        deletion = np.zeros(length, np.uint8)

        if np.sum(lower_case) > 0:
            # positions of deletions
            pos = np.where(lower_case == 1)[0]

            # shift by occurrence
            lower_case = pos - np.arange(pos.shape[0])

            # position of deletions in cleaned sequence
            # and their length
            pos, num = np.unique(lower_case, return_counts=True)

            # append to the matrix of insetions
            deletion[pos] = np.clip(num, 0, 255).astype(np.uint8)  # to save memory

        sequence = raw_sequence.translate(table)
        sequence = mapping_view.map(np.array(list(sequence)))
        sequences.append(sequence)
        deletions.append(deletion)
    query_sequence = np.array(list(query_sequence))
    sequences = np.stack(sequences).astype(np.uint8)
    deletions = np.stack(deletions).astype(np.int32)
    deletion_mean = 2 * np.arctan(deletions.astype(np.float32) / 3) / np.pi
    deletion_mean = deletion_mean.mean(axis=0).astype(np.float32)
    profile = np.eye(max_idx + 1, dtype=np.int32)[
        sequences
    ]  # for now, protein only
    profile = np.mean(profile, axis=0).astype(np.float32)

    return {
        "query_sequence": query_sequence,
        "aligned_sequences": sequences,
        "deletions": deletions,
        "deletion_mean": deletion_mean,
        "profile": profile,
    }


def parse_headers(headers: list[str]) -> dict[str, np.ndarray]:
    """Extract information from a3m FASTA headers.

    The function supports three formats:

    1. UniRef-style header:
    Example:
    >UniRef100_W5NM83 G_PROTEIN_RECEP_F1_2 domain-containing protein n=1 Tax=Lepisosteus oculatus TaxID=7918 RepID=W5NM83_LEPOC

    Extracts:
        - db_name: "UniRef100"
        - db_id:   "W5NM83"
        - species: "Lepisosteus oculatus"
        - rep_id:  "W5NM83_LEPOC"

    2. Pipe-delimited UniProt header:
    Example:
    >tr|A0A060WKI3|A0A060WKI3_ONCMY Uncharacterized protein OS=Oncorhynchus mykiss GN=GSONMT00072548001 PE=3 SV=1

    Extracts:
        - db_name: "tr"
        - db_id:   "A0A060WKI3"
        - species: "Oncorhynchus mykiss"
        - rep_id:  "A0A060WKI3_ONCMY"

    3. BFD output header:
    Example:
    >SRR4029434_2280741
    >APCry4251928276_1046603.scaffolds.fasta_scaffold646995_1 # 3 # 410 # 1 # ID=646995_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.426
    # TODO extract species info from bfd db

    Extracts:
        - db_name: "bfd"
        - db_id:   "SRR4029434_2280741"
        - species: "N/A"
        - rep_id:  "SRR4029434_2280741"

    Returns
    -------
        A dictionary with keys "db_name", "db_id", "species", and "rep_id".
    """
    # Pattern 1: UniRef-style header (with Tax=... and RepID=...)
    pattern1 = re.compile(
        r"^(?P<db_name>UniRef\d+)_"
        r"(?P<db_id>\S+).*?Tax=(?P<species>.*?)\s+TaxID=\S+\s+RepID=(?P<rep_id>\S+)",
        re.IGNORECASE,
    )

    # Pattern 2: Pipe-delimited UniProt header (with OS=...)
    pattern2 = re.compile(
        r"^(?P<db_name>[^|]+)\|"
        r"(?P<db_id>[^|]+)\|"
        r"(?P<rep_id>[^|]+)\s+.*?OS=(?P<species>.*?)\s+(?=GN=|PE=|SV=)",
        re.IGNORECASE,
    )

    result = None
    database_list = []
    database_id = []
    species_list = []
    rep_id_list = []
    for ii, header in enumerate(headers):
        if ii == 0:
            database_list.append("query")
            database_id.append("query")
            species_list.append("query")
            rep_id_list.append("query")
            continue
        for pattern in (pattern1, pattern2):
            match = pattern.search(header)
            if match:
                result = match.groupdict()
                # For pattern3, assign default values for missing keys.
                if "species" not in result or not result.get("species"):
                    result["species"] = "N/A"
                if "rep_id" not in result or not result.get("rep_id"):
                    result["rep_id"] = "N/A"
                break

        if result is not None:
            database = result.get("db_name", "N/A").lower()
            db_id = result.get("db_id", "N/A")
            species = result.get("species", "N/A")
            rep_id = result.get("rep_id", "N/A")
        else:
            # Pattern 3: BFD output header (default values).
            # Species information is not extracted from the BFD database.
            database = "bfd"
            db_id = header[1:].split()[0]  # Remove '>' and take first part
            species = "N/A"
            rep_id = db_id
        database_list.append(database)
        database_id.append(db_id)
        species_list.append(species)
        rep_id_list.append(rep_id)

    database_list = np.array(database_list, dtype="S")
    database_id = np.array(database_id, dtype="S")
    species_list = np.array(species_list, dtype="S")
    rep_id_list = np.array(rep_id_list, dtype="S")

    return {
        "database": database_list,
        "database_id": database_id,
        "species": species_list,
        "rep_id": rep_id_list,
    }


def build_dict(
    sequences: dict[str, np.ndarray],
    headers: dict[str, np.ndarray],
) -> dict[str, dict[str, np.ndarray]]:
    """Build a feature container from parsed a3m data."""
    return {
        "sequences": sequences,
        "headers": headers,
    }
