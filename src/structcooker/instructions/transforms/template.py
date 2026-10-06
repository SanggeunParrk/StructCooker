import copy
import fnmatch
import io
import os
import re
import subprocess
import time
from collections.abc import Callable
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, cast

import kalign
import numpy as np
import zstandard
from biomol.core.index import IndexTable

from structcooker.instructions.readers.io import load_bytes, load_raw_data
from structcooker.mols import CIFMol, TemplateMol
from structcooker.utils.seq_id import seq_id_from_name, seq_id_shard_path

if TYPE_CHECKING:
    from biomol.core.types import BioMolDict


def _is_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def _finished_files(root: Path, suffix: str) -> set[str]:
    """Names of the non-empty ``*<suffix>`` files anywhere under ``root``.

    One directory walk instead of a stat per expected output: on the shared filesystem,
    checking ~1.7M mostly-absent paths one by one takes tens of minutes.
    """
    found: set[str] = set()
    if not root.is_dir():
        return found
    stack = [root]
    while stack:
        with os.scandir(stack.pop()) as it:
            for entry in it:
                if entry.is_dir(follow_symlinks=False):
                    stack.append(Path(entry.path))
                elif entry.name.endswith(suffix) and entry.stat().st_size > 0:
                    found.add(entry.name)
    return found


def _run_command(
    command: list[str],
    *,
    env: dict[str, str] | None = None,
) -> None:
    result = subprocess.run(  # noqa: S603 (fixed internal tool argv)
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
    done_dir: Path | None = None,
    done_suffix: str = ".out",
) -> list[dict[str, Path]]:
    """Scan a directory recursively for A3M files matching a pattern and return a list of input/output path pairs.

    With ``done_dir``, an a3m whose final product (``<done_dir>/<shard>/<stem><done_suffix>``)
    already exists is left out of the list. Downstream steps skip finished items anyway, but
    work is striped over nodes by list position, so on a resume the unfinished items would
    stay on whichever nodes held them; dropping finished ones first spreads the rest evenly.
    """
    result = []
    done = _finished_files(Path(done_dir), done_suffix) if done_dir is not None else set()

    def _scan(dir_path: Path) -> None:
        with os.scandir(dir_path) as it:
            for entry in it:
                if entry.is_dir(follow_symlinks=False):
                    _scan(Path(entry.path))
                elif fnmatch.fnmatch(entry.name, pattern):
                    seq_id = seq_id_from_name(entry.name)
                    if f"{Path(entry.name).stem}{done_suffix}" in done:
                        continue
                    output_parent = (
                        seq_id_shard_path(output_dir, seq_id)
                        if seq_id is not None
                        else output_dir
                    )
                    result.append(
                        {
                            "input_a3m_path": Path(entry.path),
                            "output_path": output_parent
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
    print(f"HHsearch completed for {msa_path}, output saved to {hhr_path}")  # noqa: T201 (CLI progress)
    return f"Done {hhr_path.name}"


def run_hmmbuild(input_a3m_path: Path, hmm_path: Path | None) -> str:
    """Run hmmbuild to convert an a3m MSA to an HMM; return a status string."""
    # run hmmbuild to convert a3m to hmm
    if hmm_path is None:
        hmm_path = input_a3m_path.with_suffix(".hmm")
    # Written under a temporary name and renamed when complete: an interrupted run must not
    # leave a partial file that the non-empty check below would then accept as finished.
    tmp_path = hmm_path.with_name(hmm_path.name + ".tmp")
    command = [
        "hmmbuild",
        str(tmp_path),
        str(input_a3m_path),
    ]
    if hmm_path.exists() and hmm_path.stat().st_size > 0:
        print(  # noqa: T201 (CLI progress)
            f"HMM file {hmm_path} already exists and is non-empty. Skipping hmmbuild for {input_a3m_path}.",
        )
        return f"Skip {hmm_path.name} (already exists and is non-empty)"
    try:
        hmm_path.parent.mkdir(parents=True, exist_ok=True)
        _run_command(command)
        tmp_path.replace(hmm_path)
        return "hmm file created at: " + str(hmm_path)
    except Exception as e:
        msg = f"Error running hmmbuild for {input_a3m_path}: {e}"
        raise RuntimeError(msg) from e


def run_hmmsearch(
    output_dir: Path,
    hmm_path: Path,
    fasta_path: Path,
    hmmbuild_results: object = None,
) -> str:
    """Run hmmsearch of an HMM against a FASTA; return a status string.

    ``hmmbuild_results`` is unused; it only makes this step depend on run_hmmbuild so the
    executor builds the HMM (at ``hmm_path``) before searching it.
    """
    _ = hmmbuild_results
    if not _is_nonempty(hmm_path):
        msg = f"HMM file does not exist or is empty: {hmm_path}"
        raise FileNotFoundError(msg)
    if not _is_nonempty(fasta_path):
        msg = f"FASTA file does not exist or is empty: {fasta_path}"
        raise FileNotFoundError(msg)
    seq_id = seq_id_from_name(hmm_path.name)
    output_parent = (
        seq_id_shard_path(output_dir, seq_id) if seq_id is not None else output_dir
    )
    output_parent.mkdir(parents=True, exist_ok=True)
    output_path = output_parent / f"{hmm_path.stem}.out"
    if output_path.exists() and output_path.stat().st_size > 0:
        print(  # noqa: T201 (CLI progress)
            f"Output file {output_path} already exists and is non-empty. Skipping hmmsearch for {hmm_path}.",
        )
        return f"Skip {output_path.name} (already exists and is non-empty)"
    # hmmsearch writes -o as it goes; see run_hmmbuild for why it goes to a temporary name.
    tmp_path = output_path.with_name(output_path.name + ".tmp")
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
        str(tmp_path),
        str(hmm_path),
        str(fasta_path),
    ]
    try:
        _run_command(command)
        tmp_path.replace(output_path)
        return f"hmmsearch completed for {hmm_path} against {fasta_path}"  # noqa: TRY300 (return kept in try for clarity)
    except Exception as e:
        msg = f"Error running hmmsearch for {hmm_path}: {e}"
        raise RuntimeError(msg) from e


def remove_lower_from_a3m(input_a3m_path: Path, output_path: Path | None) -> str:
    """Strip lowercase (insertion) columns from an a3m; return a status string.

    Reads plain or zstd-compressed (``.zst``) a3m and drops ``#`` metadata lines (the
    ColabFold format AFDB ships). The output is written under a temporary name and renamed
    when complete, so an existing non-empty output is a finished one and is kept.
    """
    if output_path is None:
        output_path = input_a3m_path.with_suffix(".no_lower.a3m")
    if _is_nonempty(output_path):
        return f"Skip {output_path.name} (already exists and is non-empty)"
    tmp_path = output_path.with_name(output_path.name + ".tmp")
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with _open_a3m_text(input_a3m_path) as infile, tmp_path.open("w") as outfile:
            for line in infile:
                if line.startswith("#"):
                    continue
                if line.startswith(">"):
                    outfile.write(line)
                else:
                    # remove all lowercase letters from the sequence lines
                    line = "".join(c for c in line if not c.islower())  # noqa: PLW2901 (intentional in-loop rewrite)
                    outfile.write(line)
        tmp_path.replace(output_path)
        return (  # noqa: TRY300 (return kept in try for clarity)
            f"Lowercase letters removed from {input_a3m_path}, saved to {output_path}"
        )
    except Exception as e:
        msg = f"Error processing {input_a3m_path} to remove lowercase letters: {e}"
        raise RuntimeError(msg) from e


def _open_a3m_text(path: Path) -> "io.TextIOBase":
    if path.suffix != ".zst":
        return path.open("r")
    raw = path.open("rb")
    return io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(raw, closefd=True))


def list_afm_msa_wo_lower(
    msa_dir: Path,
    afm_msa_seqid_path: Path,
    output_dir: Path,
    generated_msa_dir: Path | None = None,
) -> list[dict[str, Path]]:
    """List AFDB entity MSAs to strip, one per seq_id, as <output_dir>/<shard>/<seq_id>.a3m.

    Entities sharing a sequence share a seq_id; the one kept is the smallest path, the same
    choice the msa LMDB makes (duplicate_input_policy first_path), so templates are searched
    from the MSA the DB holds. Outputs that already exist are left out, so a resume only
    lists what is left. An entity the release does not ship is read from
    ``generated_msa_dir`` (materials/intermediate), where the MSA we generated for it lives.
    """
    generated = {p.name: p for p in Path(generated_msa_dir).glob("AF-*-msa_v1.a3m.zst")} \
        if generated_msa_dir is not None else {}
    chosen: dict[str, Path] = {}
    with Path(afm_msa_seqid_path).open() as handle:
        for line in handle:
            entity, _, seq_id = line.rstrip("\n").partition("\t")
            number = entity.split("-")[1]
            name = f"{entity}-msa_v1.a3m.zst"
            path = generated.get(name) or Path(msa_dir) / number[-3:] / name
            if seq_id not in chosen or str(path) < str(chosen[seq_id]):
                chosen[seq_id] = path
    done = _finished_files(Path(output_dir), ".a3m")
    items = []
    for seq_id, path in sorted(chosen.items()):
        if f"{seq_id}.a3m" not in done:
            out = seq_id_shard_path(Path(output_dir), seq_id) / f"{seq_id}.a3m"
            items.append({"input_a3m_path": path, "output_path": out})
    return items


def parse_hmm_query_mapping(hmm_path: Path) -> tuple[dict[int, int], dict[int, int]]:
    """Parse an HMMER3 .hmm file with 'MAP yes' into query<->HMM index maps.

    Build:
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


def cif_record_to_chains(record: dict) -> dict[str, object]:
    """Split a decoded cif record into per-chain CIFMol dicts.

    Mirrors :func:`load_cifmol`'s selection for every base chain: the base chain
    is taken from the highest-total-occupancy assembly that contains it, then
    extracted. Returns ``{base_chain: biomoldict}`` so a per-chain cif LMDB can
    be prebuilt once, turning template lookups into a light keyed read.

    BioMol's ``extract()`` mutates the array buffers it works on, so each chain
    is built from a deep copy of the chosen assembly (fresh buffers) and
    extracted exactly once -- the same contract the per-hit path relies on.
    Selection reads occupancy / chain ids straight from the raw arrays to avoid
    building assemblies that are never used.
    """
    assemblies = record["assembly_dict"]
    metadata = record["metadata_dict"]
    # base chain -> (assembly total occupancy, assembly key)
    best: dict[str, tuple[float, str]] = {}
    for cif_key, item in assemblies.items():
        occupancy = np.asarray(
            item["atoms"]["nodes"]["occupancy"]["value"], dtype=np.float64,
        )
        occup_sum = float(np.nan_to_num(occupancy, nan=0.0).sum())
        chain_ids = item["chains"]["nodes"]["chain_id"]["value"]
        for base in {str(cid).split("_")[0] for cid in chain_ids}:
            if base not in best or occup_sum > best[base][0]:
                best[base] = (occup_sum, cif_key)

    chains: dict[str, object] = {}
    for base, (_, cif_key) in best.items():
        assembly_id, model_id, alt_id = cif_key.split("_")
        md = dict(metadata)
        md["assembly_id"], md["model_id"], md["alt_id"] = assembly_id, model_id, alt_id
        biomol = copy.deepcopy(assemblies[cif_key])  # fresh buffers for extract()
        biomol["metadata"] = md
        cifmol = CIFMol.from_dict(cast("BioMolDict", biomol))
        full = find_first(f"{base}_", cifmol.chains.chain_id.value)
        if full is None:
            continue
        chains[base] = cifmol.chains[cifmol.chains.chain_id == full].extract().to_dict()
    return chains


# Leave out assemblies above this atom count: a handful of giant assemblies (virus
# capsids, ribosome polysomes) would otherwise blow one worker's memory. Only the
# oversized ASSEMBLIES are left out, not the entry: a capsid's asymmetric unit is small,
# and its chains are what a template needs. (Until 2026-09-30 one oversized assembly
# dropped the whole entry -- 261 PDB entries, 47,588 chains, were missing from cif_chain.)
# env DC_CIFCHAIN_MAX_ATOMS; 0 disables it.
CIFCHAIN_MAX_ATOMS = int(os.environ.get("DC_CIFCHAIN_MAX_ATOMS", "1500000") or "0")


def _atom_count(assembly: dict) -> int:
    return len(assembly["atoms"]["nodes"]["id"]["value"])


def _within_atom_cap(assemblies: dict) -> dict:
    """Return the assemblies small enough to decode (all of them when the cap is off)."""
    if not CIFCHAIN_MAX_ATOMS:
        return assemblies
    return {k: a for k, a in assemblies.items() if _atom_count(a) <= CIFCHAIN_MAX_ATOMS}


def adapt_cif_record_to_chains(data: dict) -> dict:
    """Reader adapter for the ``chain/cif_chain`` rebuild (split_entries + explode).

    Splits a decoded cif record (``{assembly_dict, metadata_dict}``) into per-chain
    sub-entries ``{base_chain: {"biomoldict": biomoldict}}`` -- best-occupancy
    assembly per chain, via :func:`cif_record_to_chains` -- so the rebuild writes one
    record per chain keyed ``<pdbid>_<base_chain>`` (Schema C). Assemblies above
    ``CIFCHAIN_MAX_ATOMS`` are left out (the entry only if all of them are), and an
    unparsable record is skipped too.
    """
    assemblies = _within_atom_cap(data.get("assembly_dict") or {})
    if not assemblies:
        return {}
    try:
        chains = cif_record_to_chains({**data, "assembly_dict": assemblies})
    except Exception:  # noqa: BLE001 - skip unparsable records (reported as failed)
        return {}
    return {base: {"biomoldict": bd} for base, bd in chains.items()}


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

    return {
        "atoms": atom_dict,
        "residues": residue_dict,
        "chains": chain_dict,
        "index_table": index_table.to_dict(),
        "metadata": metadata,
    }


def build_chain_template(
    file_path: Path,
    template_metadata_map: dict[str, dict],
    chain2seqid: dict[str, str],
    seqid2earliest_date: dict[str, time.struct_time],
    filtered_seqid2seq: dict[str, str],
    hmm_dir: Path,
    cif_db_path: Path,
    date_cutoff: str = "2021-09-30",
    day_diff_cutoff: int = 60,
    topk: int = 20,
) -> dict:
    """Build one chain's template mols, keyed by chain instead of seq_id.

    ``file_path`` is the ``{pdbid}_{chain}`` work item. Template hits come from
    the chain's seq_id hmmsearch output but are filtered against THIS chain's own
    deposition date (not the sequence's earliest occurrence). Returns
    ``{hit: template_mol}`` (empty dict if the chain has no usable templates).
    """
    chain = Path(file_path).name
    seq_id = chain2seqid.get(chain)
    meta = template_metadata_map.get(chain)
    query_seq = filtered_seqid2seq.get(seq_id) if seq_id else None
    if seq_id is None or meta is None or query_seq is None:
        return {}
    hmm_path = Path(hmm_dir) / seq_id[0] / seq_id[-3:] / f"{seq_id}.out"
    if not hmm_path.exists():
        return {}
    template_chain_ids, _earliest, _q = extract_sequences(
        hmm_path, seqid2earliest_date, filtered_seqid2seq,
    )
    # Drop hits absent from the metadata map (chains not in the current PDB set /
    # naming mismatches); they carry no date or sequence to align against.
    template_chain_ids = [h for h in template_chain_ids if h in template_metadata_map]
    worker = filter_and_align_template_chain_ids(
        date_cutoff=date_cutoff, day_diff_cutoff=day_diff_cutoff, topk=topk,
    )
    align_results = worker(query_seq, meta["deposition_date"], template_chain_ids, template_metadata_map)
    return load_templates(Path(cif_db_path), align_results)


_HMM_HIT_RE = re.compile(
    r"""
    ^\s*
    [0-9.eE+-]+ \s+[0-9.]+ \s+[0-9.]+
    \s+[0-9.eE+-]+ \s+[0-9.]+ \s+[0-9.]+
    \s+[0-9.]+ \s+\d+
    \s+(\S+)
    \s+\|
    """,
    re.MULTILINE | re.VERBOSE,
)


def _parse_hmm_hits(text: str) -> list[str]:
    """Ordered unique template chain ids ('<pdbid>_<chain>') from hmm .out text."""
    hits: list[str] = []
    seen: set[str] = set()
    for raw in _HMM_HIT_RE.findall(text):
        tid = "_".join(raw.split("_")[:2])
        if tid not in seen:
            seen.add(tid)
            hits.append(tid)
    return hits


def build_seqid_template_mols(
    file_path: Path,
    query_seqs: dict[str, str],
    template_seqs: dict[str, str],
    reduced_hmm_dir: Path,
    cif_chain_db_path: Path,
    min_seq_len: int = 10,
    min_query_coverage: float = 0.1,
    max_query_coverage: float = 0.95,
    max_keep: int | None = None,
    allow_missing_chains: bool = False,
) -> tuple[dict, list[str]]:
    """Phase 3: build the union template mols for one seq_id.

    ``file_path`` is a seq_id. The reduced hmm (written by Phase 2) already holds
    only the union of hits selected -- with the per-chain 60-day date filter --
    by any chain of this seq_id, so here we build a mol for EVERY hit (no date
    re-filter): kalign + query-coverage + ``to_template_mol`` via the per-chain
    cif LMDB. The expensive alignment is thus done once per seq_id instead of
    once per chain.

    Only SEQUENCES are needed (not dates, not the full 16.9 M seq_id map):
    ``query_seqs`` (seq_id -> seq) and ``template_seqs`` (chain -> seq) are the
    lean maps that replace load_template_metadata's ~8 GB of state.
    With ``max_keep``, hits are aligned in e-value order only until ``max_keep`` pass
    the coverage filter, and only those are built: the final top-k after filtering, for
    a set with no per-chain step after this one (teddymer). Unset (PDB), every hit is
    built and Phase 4 selects per chain.

    ``allow_missing_chains``: a hit whose chain is absent from the per-chain LMDB is
    skipped and the next hit takes its place, instead of failing the record (PDB keeps
    the strict default). With it, hits are built in rounds until ``max_keep`` load.
    Returns ``(template_mols, template_ids)``.
    """
    seq_id = Path(file_path).name
    query_seq = query_seqs.get(seq_id)
    if not query_seq:
        return {}, []
    hmm_path = Path(reduced_hmm_dir) / seq_id[0] / seq_id[-3:] / f"{seq_id}.reduced.out"
    if not hmm_path.exists():
        return {}, []
    # Parse hits directly (extract_sequences derives the seq_id from the file
    # stem, which is '<seq_id>.reduced' here).
    hits = _parse_hmm_hits(hmm_path.read_text(encoding="utf-8"))
    mols: dict = {}
    align_results: dict[str, tuple[str, str]] = {}
    for tid in hits:
        tseq = template_seqs.get(tid)
        if tseq is None or len(tseq) < min_seq_len:
            continue
        aligned_query, aligned_template, coverage = run_kalign(
            query_seq=query_seq, template_seq=tseq,
        )
        if min_query_coverage <= coverage <= max_query_coverage:
            align_results[tid] = (aligned_query, aligned_template)
            if max_keep is not None and len(mols) + len(align_results) >= max_keep:
                if not allow_missing_chains:
                    break
                loaded, _ = load_templates_with_report(Path(cif_chain_db_path), align_results)
                mols.update(loaded)
                align_results = {}
                if len(mols) >= max_keep:
                    break
    if allow_missing_chains:
        if align_results:
            loaded, _ = load_templates_with_report(Path(cif_chain_db_path), align_results)
            mols.update(loaded)
        return mols, list(mols.keys())
    mols = load_templates_from_chain_db(Path(cif_chain_db_path), align_results)
    return mols, list(mols.keys())


def build_chain_from_seqid_db(
    file_path: Path,
    chain2seqid: dict[str, str],
    chain2templates: dict[str, list[str]],
    seqid_template_db_path: Path,
    topk: int = 20,
) -> dict:
    """Phase 4: assemble one chain's template mols by lookup (no alignment).

    ``file_path`` is a ``{pdbid}_{chain}``. ``chain2templates[chain]`` is the
    chain's date-filtered candidate ids in e-value order (Phase 2, up to ~60).
    We keep those that survived Phase 3's coverage/length filter (i.e. are in the
    seq_id union mol DB) and take the top ``topk`` in e-value order -- so coverage
    drops are backfilled from lower-ranked candidates, matching AF3's "keep up to
    20 after filtering". Pure keyed read + select -- no kalign, no CIF decode.
    """
    chain = Path(file_path).name
    seq_id = chain2seqid.get(chain)
    tids = chain2templates.get(chain)
    if seq_id is None or not tids:
        return {}
    raw = load_raw_data(seq_id, cast("Path", str(seqid_template_db_path)))
    if raw is None:
        return {}
    union = load_bytes(raw).get("template_mols", {})
    kept = [tid for tid in tids if tid in union][:topk]
    return {tid: union[tid] for tid in kept}


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
        try:
            cifmol = load_cifmol(cif_db_path, pdb_id.lower(), chain_id)
            chain_id = find_first(f"{chain_id}_", cifmol.chains.chain_id.value)
            cifmol = cifmol.chains[cifmol.chains.chain_id == chain_id].extract()
            template_mols[full_id] = to_template_mol(cifmol, align_result)
        except Exception:  # noqa: BLE001,S112 - skip templates missing from the CIF DB
            continue
    return template_mols


def rank_template_hits_by_coverage(
    template_hits: dict[str, dict[str, object]],
    query_len: int,
    min_coverage: float = 0.1,
) -> dict[str, dict[str, object]]:
    """Order template hits by query coverage (descending), dropping low-coverage ones.

    Coverage = unique query residues aligned (``idx_map[:, 0]``) / ``query_len``.
    Computed from ``idx_map`` alone -- no per-chain decode -- so the caller can
    decode the highest-coverage hits first and stop once it has enough
    (see :func:`load_templates_from_chain_db`'s ``max_keep``).
    """
    if not template_hits or query_len <= 0:
        return template_hits
    scored: list[tuple[float, str]] = []
    for hit, payload in template_hits.items():
        idx_map = np.asarray(payload["idx_map"], dtype=np.int64)
        coverage = len(np.unique(idx_map[:, 0])) / query_len
        if coverage >= min_coverage:
            scored.append((min(coverage, 1.0), hit))
    scored.sort(key=lambda x: x[0], reverse=True)
    return {hit: template_hits[hit] for _, hit in scored}


def load_templates_with_report(
    cif_chain_db_path: Path,
    align_results: dict[str, tuple[str, str]],
    max_keep: int | None = None,
) -> tuple[dict, dict]:
    """Load ranked hits, recording absent and misaligned chains; corrupt hits fail the record.

    A hit whose alignment points past the end of our copy of the chain (the aligner saw a
    longer chain than the one in cif_chain, e.g. 8t0v_C: residue index 259 of 259) is
    skipped and recorded in ``misaligned_hits``, and the next hit takes its place -- one
    stale hit must not cost the record all its templates. A chain that fails to decode is
    still an error.
    """
    template_mols: dict = {}
    missing: list[str] = []
    misaligned: list[str] = []
    unselected: list[str] = []
    for full_id, align_result in align_results.items():
        if max_keep is not None and len(template_mols) >= max_keep:
            unselected.append(full_id)
            continue
        pdb_id, chain_id = full_id.split("_", 1)
        raw = load_raw_data(f"{pdb_id.lower()}_{chain_id}", cif_chain_db_path)
        if raw is None:
            missing.append(full_id)
            continue
        try:
            cifmol = CIFMol.from_dict(cast("BioMolDict", load_bytes(raw)))
        except Exception as exc:
            msg = f"Template hit {full_id} failed decoding"
            raise ValueError(msg) from exc
        try:
            template_mols[full_id] = to_template_mol(cifmol, align_result)
        except IndexError:
            misaligned.append(full_id)
    report = {"candidate_count": len(align_results), "loaded_hits": list(template_mols),
              "missing_chain_hits": missing, "misaligned_hits": misaligned,
              "not_selected_hits": unselected}
    return template_mols, report


def load_templates_from_chain_db(
    cif_chain_db_path: Path,
    align_results: dict[str, tuple[str, str]],
    max_keep: int | None = None,
) -> dict:
    """Strict compatibility API; use the reported API to allow absent reference chains."""
    mols, report = load_templates_with_report(cif_chain_db_path, align_results, max_keep)
    if report["missing_chain_hits"]:
        msg = f"Missing reference template chains: {report['missing_chain_hits']}"
        raise ValueError(msg)
    return mols


# Unittest functions


def length_check(
    query_seq_id: str,
    seqid2seq: dict[str, list[str]],
    templatemol_dict: dict[str, TemplateMol],
) -> str | None:
    """Unittest function to check if the length of the query sequence matches the length of the template mol sequences."""
    query_seq = seqid2seq.get(query_seq_id)
    if not query_seq:
        return f"Query sequence for ID {query_seq_id} not found."
    query_seq = query_seq[0] if isinstance(query_seq, list) else query_seq
    for full_id, template_mol in templatemol_dict.items():
        template_seq_len = len(template_mol.residues)
        if len(query_seq) != template_seq_len:
            return f"Length mismatch for {full_id}: query length {len(query_seq)} != template length {template_seq_len}"
    return None


def adapt_cif_record_to_chain_inputs(data: dict) -> dict:
    """Prepare chain inputs once, with deterministic maximum-occupancy candidates.

    Candidates come only from assemblies within ``CIFCHAIN_MAX_ATOMS``; an entry is
    skipped only when every assembly is above it.
    """
    assemblies = _within_atom_cap(data.get("assembly_dict") or {})
    if not assemblies:
        return {}
    best: dict[str, tuple[float, str]] = {}
    # Numeric assembly/model order makes ties independent of dict/hash order.
    ordered = sorted(assemblies, key=lambda key: tuple(
        (0, int(part)) if part.isdecimal() else (1, part)
        for part in key.split("_")
    ))
    for key in ordered:
        item = assemblies[key]
        score = float(np.nan_to_num(np.asarray(
            item["atoms"]["nodes"]["occupancy"]["value"], dtype=float,
        ), nan=0.0).sum())
        for chain in item["chains"]["nodes"]["chain_id"]["value"]:
            base = str(chain).split("_")[0]
            if base not in best or score > best[base][0]:
                best[base] = score, key
    model_cache: dict = {}
    return {chain: {"record": data, "chain_id": chain, "cif_key": key,
                    "model_cache": model_cache}
            for chain, (_, key) in best.items()}


def extract_selected_chain(
    record: dict,
    chain_id: str,
    cif_key: str,
    reference_selection_path: Path | None = None,
    model_cache: dict | None = None,
) -> dict:
    """Extract a chain, optionally preserving a reference dataset's model choice."""
    if reference_selection_path is not None and _selection_reference_exists(reference_selection_path):
        from datacooker.lmdb import read_lmdb_raw

        pdbid = str(record["metadata_dict"]["id"][0]).lower()
        selected = read_lmdb_raw(reference_selection_path, f"{pdbid}_{chain_id}")
        # A reference choice pointing at an assembly above the atom cap is not honoured:
        # decoding it is what the cap exists to prevent.
        if selected is not None and selected.decode() in _within_atom_cap(record["assembly_dict"]):
            cif_key = selected.decode()
    cache_state = model_cache if model_cache is not None else {}
    if cache_state.get("record") is not record or cache_state.get("cif_key") != cif_key:
        assembly_id, model_id, alt_id = cif_key.split("_")
        biomol = copy.deepcopy(record["assembly_dict"][cif_key])
        biomol["metadata"] = {**record["metadata_dict"], "assembly_id": assembly_id,
                              "model_id": model_id, "alt_id": alt_id}
        # This cache belongs to one adapter result, so it dies with the record.
        # Retain only the most recent model, even when references choose many models.
        cache_state.clear()
        cache_state.update(record=record, cif_key=cif_key,
                           cifmol=CIFMol.from_dict(cast("BioMolDict", biomol)))
    cifmol = cache_state["cifmol"]
    full = find_first(f"{chain_id}_", cifmol.chains.chain_id.value)
    if full is None:
        msg = f"Selected model {cif_key} does not contain chain {chain_id}."
        raise ValueError(msg)
    return cifmol.chains[cifmol.chains.chain_id == full].extract().to_dict()


@cache
def _selection_reference_exists(path: Path) -> bool:
    return path.exists()
