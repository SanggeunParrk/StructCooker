"""Small fasta/TSV -> list projections feeding the template pipeline.

Tiny side-effect writers (like the precompute steps): each reads a fasta/TSV and
writes a derived list/fasta, returning a status string. Driven by materialize
configs whose ``output_data_path`` is the build-done marker.
"""
from __future__ import annotations

from pathlib import Path

_L_TYPE = "polypeptide(L)"


def _chain_id(header: str) -> str:
    """``100D_A_. | polypeptide(L) | ...`` -> ``100D_A`` (pdbid_chain)."""
    tok = header.split("|", 1)[0].strip()
    parts = tok.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else tok


def filter_fasta_polypeptide_l(fasta_path: str | Path, out_path: str | Path) -> str:
    """Write only the polypeptide(L) records of ``fasta_path`` (the hmmsearch template DB)."""
    n = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    keep = False
    with Path(fasta_path).open(encoding="utf-8") as f, Path(out_path).open("w", encoding="utf-8") as out:
        for line in f:
            if line.startswith(">"):
                keep = _L_TYPE in line
                if keep:
                    out.write(line)
                    n += 1
            elif keep:
                out.write(line)
    return f"polypeptide(L) fasta: {n} records -> {out_path}"


def fasta_chain_list(fasta_path: str | Path, out_path: str | Path) -> str:
    """Write the pdbid_chain id of every polypeptide(L) record, one per line.

    This is the protein-chain work list (``template_chain_filelist``) that the
    Phase 1/2 candidate precompute and the Phase 4 keyed build iterate over.
    """
    n = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(fasta_path).open(encoding="utf-8") as f, Path(out_path).open("w", encoding="utf-8") as out:
        for line in f:
            if line.startswith(">") and _L_TYPE in line:
                out.write(_chain_id(line[1:]) + "\n")
                n += 1
    return f"chain list: {n} polypeptide(L) chains -> {out_path}"


def tsv_key_list(tsv_path: str | Path, out_path: str | Path) -> str:
    """Write column-0 keys of a TSV, one per line (e.g. seqid_to_chains -> seq_ids)."""
    n = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(tsv_path).open(encoding="utf-8") as f, Path(out_path).open("w", encoding="utf-8") as out:
        for line in f:
            key = line.split("\t", 1)[0].strip()
            if key:
                out.write(key + "\n")
                n += 1
    return f"key list: {n} keys -> {out_path}"


def glob_file_list(sources: list[dict], out_path: str | Path) -> str:
    """Write every file matching ``{dir, glob[, recursive]}`` across ``sources``, sorted, one per line.

    For a DB whose inputs live in more than one tree -- a release's files under
    ``materials/raw`` plus files we generated under ``materials/intermediate`` -- so the
    raw tree never has to hold anything that was not downloaded.
    """
    paths: list[str] = []
    for source in sources:
        root, pattern = Path(source["dir"]), str(source.get("glob", "*"))
        found = root.rglob(pattern) if source.get("recursive") else root.glob(pattern)
        paths.extend(str(p) for p in found if p.is_file())
    paths.sort()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text("".join(f"{p}\n" for p in paths), encoding="utf-8")
    return f"file list: {len(paths)} files from {len(sources)} sources -> {out_path}"

