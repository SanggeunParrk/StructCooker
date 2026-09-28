import os
import re
import shutil
import string
import subprocess
from pathlib import Path
from typing import TypeVar

import numpy as np

from structcooker.utils.mapping import ResidueMapping
from structcooker.utils.seq_id import seq_id_shard_path

DEFAULT_DB_UR30 = Path(
    "/data/shared/cssb_data/db_protSeq/uniref30/2022_02/UniRef30_2022_02",
)
DEFAULT_DB_BFD = Path(
    "/data/shared/cssb_data/db_protSeq/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
)
DEFAULT_HHSUITE_BIN = Path("/software/hhsuite/build/bin")

# RNA MSA (AlphaFold3-style): nhmmer over Rfam / RNAcentral / nucleotide
# collection, hits realigned to the query with hmmalign, cropped to 5000.
DEFAULT_DB_RFAM = Path("/data/psk6950/rmsa_db/rfam.fasta")
DEFAULT_DB_RNACENTRAL = Path("/data/psk6950/rmsa_db/rnacentral.fasta")
DEFAULT_DB_NT = Path("/data/psk6950/rmsa_db/nucleotide_collection.fasta")
RNA_MSA_MAX_SEQUENCES = 5000
_SHORT_RNA_LEN = 50
# Long rRNA are the pathological case: scan cost scales with query length, and a
# 2-3knt rRNA against the 27GB DB set takes hours. For such queries we skip the
# 10.7GB nucleotide-collection DB (Rfam + RNAcentral already cover conserved
# rRNA well). We only do this for queries that both look like rRNA (their
# RNAcentral hits are majority rRNA-typed) AND are long enough to be expensive.
_NT_SKIP_MIN_LEN = 1000
_RRNA_HIT_FRACTION = 0.5
# Long RNAs (>= this) in the PDB are essentially all rRNA. Even after dropping
# nt, scanning the 16GB RNAcentral DB with a ~1-5knt query is O(query x DB) and
# takes many hours (a 1.5knt query at a few cores can run ~half a day), which
# stalls the whole run. For these we search Rfam only (215MB, curated rRNA
# families) -- minutes per query -- which is ample coverage for conserved rRNA.
_RFAM_ONLY_MIN_LEN = 1000
# hmmalign is single-threaded and roughly O(n_seqs x query_len^2); for very long
# rRNA (~2-3knt) with the full 5000 hits it runs for hours or effectively hangs.
# Cap the alignment depth for these -- 1000 sequences is still a deep rRNA MSA.
_ALIGN_LONG_LEN = 2000
_ALIGN_LONG_MAX_SEQS = 1000
_RNA_MSA_ALLOWED_CHARS = set("ACGUNX")


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
    result = subprocess.run(command, capture_output=True, text=True, env=env, check=False)  # noqa: S603 (fixed internal tool argv)
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
    """Write ``sequence`` to a per-seqid FASTA; return (fasta_path, out_dir)."""
    out_dir = seq_id_shard_path(output_dir, seqid) / seqid
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = out_dir / f"{seqid}.fasta"
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with fasta_path.open("w") as f:
        f.write(f">{seqid}\n{sequence}\n")
    return fasta_path, out_dir

def run_signalp(
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

def run_msa_search(
    input_fasta: Path,
    out_dir: Path,
    *,
    cpu: int = 4,
    mem: int = 20,
    db_ur30: Path,
    db_bfd: Path,
    hhsuite_bin_dir: Path,
) -> str:
    """Run the HHblits MSA search pipeline for one FASTA; return a status string."""
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

            if n75 > 2000:
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(id90cov75_file, msa0_file)
                    break
            elif n50 > 4000:  # noqa: SIM102 (branch selection differs if collapsed)
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
            if n75 > 2000:
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(bfd_id90cov75_file, msa0_file)
            elif n50 > 4000:  # noqa: SIM102 (branch selection differs if collapsed)
                if not _is_nonempty(msa0_file):
                    shutil.copyfile(bfd_id90cov50_file, msa0_file)

        if not _is_nonempty(msa0_file):
            shutil.copyfile(prev_a3m, msa0_file)

    return f"Done {input_fasta}"


def _nhmmer_command(
    query_fasta: Path,
    db_path: Path,
    tbl_out: Path,
    *,
    cpu: int,
    f3: str,
) -> list[str]:
    """AlphaFold3 RNA nhmmer invocation (E=1e-3, watson strand, RNA alphabet).

    Only the compact ``--tblout`` hit table is produced (one line per hit,
    already sorted best-E-value first) -- crucially NOT ``-A``, whose full
    per-hit alignment output is unbounded and explodes to hundreds of GB for
    conserved RNAs (e.g. rRNA) that hit millions of database sequences.
    """
    return [
        "nhmmer",
        "-o",
        "/dev/null",
        "--noali",
        "--tblout",
        str(tbl_out),
        "-E",
        "0.001",
        "--incE",
        "0.001",
        "--rna",
        "--watson",
        "--F3",
        f3,
        "--cpu",
        str(cpu),
        str(query_fasta),
        str(db_path),
    ]


def _read_top_tblout(tbl_path: Path, limit: int) -> list[tuple[float, str, str, str]]:
    """Read the top ``limit`` hits (E-value, target, from, to) from a tblout.

    nhmmer writes tblout sorted by ascending E-value, so the first ``limit``
    data rows are the best hits and the rest of the (potentially multi-GB) file
    can be ignored.
    """
    hits: list[tuple[float, str, str, str]] = []
    with tbl_path.open("r") as f:
        for raw_line in f:
            if raw_line.startswith("#"):
                continue
            cols = raw_line.split()
            if len(cols) < 13:
                continue
            target, alifrom, alito, evalue = cols[0], cols[6], cols[7], cols[12]
            try:
                hits.append((float(evalue), target, alifrom, alito))
            except ValueError:
                continue
            if len(hits) >= limit:
                break
    return hits


def _a2m_to_a3m_line(a2m_seq: str) -> str:
    """Convert one a2m aligned sequence to a3m (drop insert-state gap dots)."""
    return a2m_seq.replace(".", "")


def _parse_fasta_alignment(path: Path) -> dict[str, str]:
    """Parse an aligned-FASTA / a2m file into {name: aligned sequence}."""
    aligned: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with path.open("r") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    aligned[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        aligned[name] = "".join(chunks)
    return aligned


def _nhmmer_top_hits(
    input_fasta: Path,
    work: Path,
    dbs: tuple[tuple[str, Path], ...],
    *,
    cpu: int,
    f3: str,
    max_sequences: int,
) -> dict[str, list[tuple[float, str, str, str]]]:
    """Run nhmmer over each DB; return best hits grouped by DB tag.

    Each value is a list of ``(evalue, target, from, to)`` for that DB's best
    hits (already E-value ordered), capped at ``max_sequences`` per DB.
    """
    hits_by_db: dict[str, list[tuple[float, str, str, str]]] = {}
    for tag, db_path in dbs:
        tbl_out = work / f"{tag}.tbl"
        if not _is_nonempty(tbl_out):
            _run_command(_nhmmer_command(input_fasta, db_path, tbl_out, cpu=cpu, f3=f3))
        hits_by_db[tag] = _read_top_tblout(tbl_out, max_sequences)
    return hits_by_db


def _esl_sfetch_subseqs(
    db_path: Path,
    hits: list[tuple[float, str, str, str]],
    namefile: Path,
    out_fasta: Path,
) -> dict[str, str]:
    """Fetch hit subsequences from ``db_path`` via esl-sfetch; return {name: seq}.

    ``namefile`` is written in Easel GDF format (``newname from to source``)
    and the DB is expected to carry an ``.ssi`` index built once up front with
    ``esl-sfetch --index``.
    """
    if not hits:
        return {}
    with namefile.open("w") as f:
        for _evalue, target, a, b in hits:
            f.write(f"{target}/{a}-{b} {a} {b} {target}\n")
    _run_command(
        ["esl-sfetch", "-Cf", "-o", str(out_fasta), str(db_path), str(namefile)],
    )
    return _read_fasta(out_fasta)


def _read_fasta(path: Path) -> dict[str, str]:
    """Parse a (unaligned) FASTA file into {name: sequence}."""
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with path.open("r") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def _search_and_fetch(
    input_fasta: Path,
    work: Path,
    dbs: tuple[tuple[str, Path], ...],
    *,
    cpu: int,
    f3: str,
    max_sequences: int,
    raw_hits: dict[str, str],
    evalues: dict[str, float],
) -> None:
    """Run nhmmer over ``dbs``, fetch hit subsequences, update raw_hits/evalues."""
    db_paths = dict(dbs)
    hits_by_db = _nhmmer_top_hits(
        input_fasta, work, dbs, cpu=cpu, f3=f3, max_sequences=max_sequences,
    )
    for tag, hits in hits_by_db.items():
        for evalue, target, a, b in hits:
            evalues.setdefault(f"{target}/{a}-{b}", evalue)
        fetched = _esl_sfetch_subseqs(
            db_paths[tag], hits, work / f"{tag}.gdf", work / f"{tag}.hits.fasta",
        )
        for name, seq in fetched.items():
            raw_hits.setdefault(name, seq.upper())


def _rrna_fraction(hits_fasta: Path) -> tuple[int, int]:
    """Return (rRNA-typed hit count, total hit count) for a fetched hits FASTA.

    RNAcentral FASTA headers carry the RNA type as the token after the name
    (e.g. ``>URS0000AF1685/12-83 rRNA from 3 species``); esl-sfetch preserves it,
    so a query whose hits are majority ``rRNA`` is itself rRNA.
    """
    rrna = total = 0
    if not hits_fasta.exists():
        return 0, 0
    with hits_fasta.open("r") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            total += 1
            parts = line.split()
            if len(parts) > 1 and parts[1] == "rRNA":
                rrna += 1
    return rrna, total


def _read_query_sequence(input_fasta: Path) -> str:
    query_seq = ""
    with input_fasta.open("r") as f:
        for line in f:
            if not line.startswith(">"):
                query_seq += line.strip()
    return query_seq


def _sanitize_rna_query_sequence(seq: str) -> str:
    """Map non-RNA polymer letters, e.g. aminoacylated PHE/P, to X."""
    return "".join(nt if nt in _RNA_MSA_ALLOWED_CHARS else "X" for nt in seq.upper())


def run_rna_msa_search(
    input_fasta: Path,
    out_dir: Path,
    *,
    cpu: int = 8,
    db_rfam: Path = DEFAULT_DB_RFAM,
    db_rnacentral: Path = DEFAULT_DB_RNACENTRAL,
    db_nt: Path = DEFAULT_DB_NT,
    max_sequences: int = RNA_MSA_MAX_SEQUENCES,
    nt_skip_min_len: int = _NT_SKIP_MIN_LEN,
    rfam_only_min_len: int = _RFAM_ONLY_MIN_LEN,
) -> str:
    """Build an AlphaFold3-style RNA MSA for one query FASTA.

    nhmmer searches Rfam and RNAcentral (and, unless the query is a long rRNA,
    the nucleotide collection); the pooled hits are fetched with esl-sfetch,
    realigned to the query with hmmalign, ranked by E-value, deduplicated by
    ``accession/from-to`` and cropped to ``max_sequences``. The result is
    written as ``<seqid>.a3m`` in ``out_dir`` (query first). If no hits are
    found a query-only a3m is written.

    Long rRNA (query >= ``nt_skip_min_len`` whose RNAcentral hits are majority
    rRNA-typed) skip the huge nucleotide-collection scan, which would otherwise
    take hours; every other query -- including long non-rRNA -- keeps the full
    three-database AlphaFold3 search.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    seqid = input_fasta.stem
    final_a3m = out_dir / f"{seqid}.a3m"
    if _is_nonempty(final_a3m):
        return f"Done {input_fasta} (cached)"

    raw_query_seq = _read_query_sequence(input_fasta)
    query_seq = _sanitize_rna_query_sequence(raw_query_seq)
    if query_seq != raw_query_seq:
        with input_fasta.open("w") as f:
            f.write(f">{seqid}\n{query_seq}\n")
    f3 = "0.02" if len(query_seq) < _SHORT_RNA_LEN else "0.00005"

    work = out_dir / "nhmmer"
    work.mkdir(parents=True, exist_ok=True)

    raw_hits: dict[str, str] = {}
    evalues: dict[str, float] = {}

    # Very long RNA (large rRNA): Rfam only -- RNAcentral/nt scans would take
    # hours and these are curated rRNA families in Rfam anyway.
    if len(query_seq) >= rfam_only_min_len:
        _search_and_fetch(
            input_fasta, work, (("rfam", db_rfam),),
            cpu=cpu, f3=f3, max_sequences=max_sequences,
            raw_hits=raw_hits, evalues=evalues,
        )
    else:
        # Stage 1: Rfam + RNAcentral. RNAcentral hit types tell us whether the
        # query is rRNA, which decides whether to run the expensive nt scan.
        _search_and_fetch(
            input_fasta, work, (("rfam", db_rfam), ("rnacentral", db_rnacentral)),
            cpu=cpu, f3=f3, max_sequences=max_sequences,
            raw_hits=raw_hits, evalues=evalues,
        )
        rrna, total = _rrna_fraction(work / "rnacentral.hits.fasta")
        is_long_rrna = (
            len(query_seq) >= nt_skip_min_len
            and total > 0
            and rrna / total >= _RRNA_HIT_FRACTION
        )
        # Stage 2: nucleotide collection, skipped only for long rRNA.
        if not is_long_rrna:
            _search_and_fetch(
                input_fasta, work, (("nt", db_nt),),
                cpu=cpu, f3=f3, max_sequences=max_sequences,
                raw_hits=raw_hits, evalues=evalues,
            )

    if not raw_hits:
        with final_a3m.open("w") as f:
            f.write(f">{seqid}\n{query_seq}\n")
        return f"Done {input_fasta} (query-only)"

    align_cap = (
        _ALIGN_LONG_MAX_SEQS if len(query_seq) >= _ALIGN_LONG_LEN else max_sequences
    )
    ordered = sorted(raw_hits, key=lambda n: evalues.get(n, float("inf")))[:align_cap]
    hits_fasta = work / "hits.fasta"
    with hits_fasta.open("w") as f:
        for name in ordered:
            f.write(f">{name}\n{raw_hits[name]}\n")

    query_hmm = work / "query.hmm"
    _run_command(["hmmbuild", "--rna", str(query_hmm), str(input_fasta)])
    aligned_a2m = work / "aligned.a2m"
    _run_command(
        [
            "hmmalign", "--rna", "--trim", "--outformat", "a2m",
            "-o", str(aligned_a2m), str(query_hmm), str(hits_fasta),
        ],
    )
    aligned = _parse_fasta_alignment(aligned_a2m)

    with final_a3m.open("w") as f:
        f.write(f">{seqid}\n{query_seq}\n")
        for hit_name in ordered:
            a2m_seq = aligned.get(hit_name)
            if a2m_seq:
                f.write(f">{hit_name}\n{_a2m_to_a3m_line(a2m_seq)}\n")
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
