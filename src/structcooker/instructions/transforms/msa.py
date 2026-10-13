import re
import string
import subprocess
from pathlib import Path
from typing import TypeVar

import numpy as np

from structcooker.utils.mapping import ResidueMapping

# Long rRNA are the pathological case: scan cost scales with query length, and a
# 2-3knt rRNA against the 27GB DB set takes hours. For such queries we skip the
# 10.7GB nucleotide-collection DB (Rfam + RNAcentral already cover conserved
# rRNA well). We only do this for queries that both look like rRNA (their
# RNAcentral hits are majority rRNA-typed) AND are long enough to be expensive.
# Long RNAs (>= this) in the PDB are essentially all rRNA. Even after dropping
# nt, scanning the 16GB RNAcentral DB with a ~1-5knt query is O(query x DB) and
# takes many hours (a 1.5knt query at a few cores can run ~half a day), which
# stalls the whole run. For these we search Rfam only (215MB, curated rRNA
# families) -- minutes per query -- which is ample coverage for conserved rRNA.
# hmmalign is single-threaded and roughly O(n_seqs x query_len^2); for very long
# rRNA (~2-3knt) with the full 5000 hits it runs for hours or effectively hangs.
# Cap the alignment depth for these -- 1000 sequences is still a deep rRNA MSA.


def _is_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def _run_command(
    command: list[str],
    *,
    env: dict[str, str] | None = None,
) -> None:
    result = subprocess.run(command, capture_output=True, text=True, env=env, check=False)  # noqa: S603 (fixed internal tool argv)
    if result.returncode != 0:
        cmd = " ".join(command)
        msg = (
            f"Command failed ({result.returncode}): {cmd}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        raise RuntimeError(msg)


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

    if not raw_sequences:
        msg = "MSA must contain a query sequence."
        raise ValueError(msg)
    query_sequence = raw_sequences[0]
    length = len(query_sequence.translate(table))
    if not length:
        msg = "MSA query has no aligned columns."
        raise ValueError(msg)
    sequences = np.empty((len(raw_sequences), length), dtype=np.uint8)
    deletions = np.zeros((len(raw_sequences), length), dtype=np.int32)
    for row, raw_sequence in enumerate(raw_sequences):
        sequence = raw_sequence.translate(table)
        if len(sequence) != length:
            msg = f"MSA row {row} has {len(sequence)} aligned columns; expected {length}."
            raise ValueError(msg)
        if any(not ("A" <= c <= "Z" or c == "-") for c in sequence):
            msg = f"MSA row {row} contains an invalid aligned character."
            raise ValueError(msg)
        sequences[row] = mapping_view.map(np.array(list(sequence)))
        positions = np.fromiter(
            (i for i, c in enumerate(raw_sequence) if "a" <= c <= "z"),
            dtype=np.int64,
        )
        if positions.size:
            columns, counts = np.unique(positions - np.arange(positions.size), return_counts=True)
            # A3M insertions belong to the following aligned column. Trailing
            # insertions have no following column and do not enter the matrix.
            keep = columns < length
            deletions[row, columns[keep]] = np.minimum(counts[keep], 255)
    deletion_mean, profile = msa_statistics(sequences, deletions, max_idx + 1)
    return {
        "query_sequence": np.array(list(query_sequence)),
        "aligned_sequences": sequences,
        "deletions": deletions,
        "deletion_mean": deletion_mean,
        "profile": profile,
    }


def msa_statistics(
    sequences: np.ndarray,
    deletions: np.ndarray,
    n_classes: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute MSA statistics without a depth x length x alphabet temporary."""
    n_rows, length = sequences.shape
    if not n_rows or not length or deletions.shape != sequences.shape:
        msg = "MSA and deletion matrices must have matching, nonempty dimensions."
        raise ValueError(msg)
    deletion_mean = (2 * np.arctan(deletions.astype(np.float32) / 3) / np.pi).mean(axis=0).astype(np.float32)
    counts = np.zeros(length * n_classes, dtype=np.int64)
    offsets = np.arange(length) * n_classes
    rows_per_chunk = max(1, 1_000_000 // length)
    for start in range(0, n_rows, rows_per_chunk):
        indices = sequences[start:start + rows_per_chunk].astype(np.int64) + offsets
        counts += np.bincount(indices.ravel(), minlength=counts.size)
    profile = (counts.reshape(length, n_classes) / n_rows).astype(np.float32)
    return deletion_mean, profile


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

    database_list = []
    database_id = []
    species_list = []
    rep_id_list = []
    for ii, raw_header in enumerate(headers):
        result = None
        header = raw_header.removeprefix(">").strip()
        if not header:
            msg = f"Empty MSA header at row {ii}."
            raise ValueError(msg)
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
            db_id = header.split()[0]
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
