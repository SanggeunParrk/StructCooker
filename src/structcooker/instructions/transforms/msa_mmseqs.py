"""MMseqs2 MSA search: one batched search per chunk of queries, ColabFold-style.

The HHblits path (``transforms/msa.py``) searches one sequence at a time, which is the
only way HHblits works. MMseqs2 must not be driven that way -- its cost is dominated by
loading the target database, so a per-sequence loop would be slower than HHblits, not
faster. Instead each work item is a CHUNK of queries that share one search.

Target is uniref30_2302 alone, no metagenomic database. Measured on 716 Teddymer domains
already run through the HHblits cascade (2026-09-23): uniref30 alone cleared the depth
threshold for 83.7%, and the resulting MSAs have a median depth of 9,925 -- far past the
2,048 cap training applies. BFD only deepens alignments that are already truncated. This
also matches how the AFDB complex set was built (MMseqs2-GPU, UniRef30 2302, --use-env 0).

Stages are the ColabFold search: search -> expandaln -> align -> result2msa -> unpackdb.
``expandaln`` is what makes a UniRef30 hit stand for its whole cluster, so the depth comes
out comparable to searching the unclustered database.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from structcooker.utils.seq_id import seq_id_shard_path

# Search sensitivity/limits. These follow ColabFold's uniref search rather than MMseqs2
# defaults: three profile iterations and a wide max-seqs are what give an a3m deep enough
# to be worth capping at 2,048.
_NUM_ITERATIONS = 3
_SENSITIVITY = 8
_MAX_SEQS = 10000
_EVALUE = 0.1


def _run(command: list[str]) -> None:
    result = subprocess.run(command, capture_output=True, text=True, check=False)  # noqa: S603
    if result.returncode != 0:
        msg = f"{command[0]} {command[1]} failed ({result.returncode}): {result.stderr[-2000:]}"
        raise RuntimeError(msg)


def write_chunk_fasta(chunk: list[tuple[str, str]], work_dir: Path) -> Path:
    """Write one chunk's queries to a FASTA keyed by seq_id."""
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta = work_dir / "query.fasta"
    with fasta.open("w") as handle:
        for seqid, sequence in chunk:
            handle.write(f">{seqid}\n{sequence}\n")
    return fasta


def run_mmseqs_msa_search(
    chunk: list[tuple[str, str]],
    chunk_index: int,
    work_dir: Path,
    output_dir: Path,
    *,
    db_uniref: Path,
    mmseqs_bin: Path,
    threads: int = 8,
) -> dict[str, int]:
    """Search one chunk and unpack an a3m per seq_id. Returns counts."""
    # Per chunk, because the scratch is wiped on entry: a shared directory would let
    # concurrent chunks delete each other's databases mid-search.
    work_dir = Path(work_dir) / f"chunk_{chunk_index:05d}"
    # Only queries whose a3m is not already on disk; the whole pipeline must stay
    # resumable, and a chunk is re-run whole if it was interrupted.
    pending = [
        (seqid, seq)
        for seqid, seq in chunk
        if not _a3m_path(output_dir, seqid).exists()
    ]
    if not pending:
        return {"queries": len(chunk), "searched": 0, "written": 0}

    if work_dir.exists():
        shutil.rmtree(work_dir)  # a partial previous attempt must not be reused
    fasta = write_chunk_fasta(pending, work_dir)
    mm = str(mmseqs_bin)
    qdb, res, exp, aln, msa = (work_dir / n for n in ("qdb", "res", "exp", "aln", "msa"))
    tmp = work_dir / "tmp"
    common = ["--threads", str(threads)]

    _run([mm, "createdb", str(fasta), str(qdb)])
    _run([
        mm, "search", str(qdb), str(db_uniref), str(res), str(tmp),
        "--num-iterations", str(_NUM_ITERATIONS), "-s", str(_SENSITIVITY),
        "-e", str(_EVALUE), "--max-seqs", str(_MAX_SEQS), "-a",
        "--db-load-mode", "2", *common,
    ])
    # A UniRef30 hit represents a cluster; expandaln pulls in the cluster's members so the
    # alignment depth matches an unclustered search.
    _run([
        mm, "expandaln", str(qdb), f"{db_uniref}_seq", str(res), f"{db_uniref}_aln",
        str(exp), "--expansion-mode", "0", "-e", "inf", "--db-load-mode", "2", *common,
    ])
    _run([
        mm, "align", str(qdb), f"{db_uniref}_seq", str(exp), str(aln),
        "-e", "10", "--max-accept", str(_MAX_SEQS), "--alt-ali", "10", "-a",
        "--db-load-mode", "2", *common,
    ])
    _run([
        mm, "result2msa", str(qdb), f"{db_uniref}_seq", str(aln), str(msa),
        "--msa-format-mode", "6", "--db-load-mode", "2", *common,
    ])
    unpack = work_dir / "unpacked"
    _run([mm, "unpackdb", str(msa), str(unpack), "--unpack-suffix", ".a3m"])

    # unpackdb names files by the database key (0.a3m, 1.a3m, ...), not by the query.
    # qdb.lookup is the key -> header mapping, and the header is the seq_id.
    key_to_seqid = {}
    with (work_dir / "qdb.lookup").open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                key_to_seqid[parts[0]] = parts[1]

    written = 0
    for src in unpack.glob("*.a3m"):
        seqid = key_to_seqid.get(src.stem)
        if seqid is None:
            msg = f"{src.name} has no entry in qdb.lookup"
            raise RuntimeError(msg)
        dst = _a3m_path(output_dir, seqid)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))
        written += 1
    shutil.rmtree(work_dir, ignore_errors=True)
    return {"queries": len(chunk), "searched": len(pending), "written": written}


def _a3m_path(output_dir: Path, seqid: str) -> Path:
    """Return the a3m path, in the layout the HHblits path and db/msa/a3m.yaml use."""
    return seq_id_shard_path(Path(output_dir), seqid) / seqid / f"{seqid}.a3m"
