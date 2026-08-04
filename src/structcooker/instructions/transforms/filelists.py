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
    with Path(fasta_path).open() as f, Path(out_path).open("w") as out:
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
    with Path(fasta_path).open() as f, Path(out_path).open("w") as out:
        for line in f:
            if line.startswith(">") and _L_TYPE in line:
                out.write(_chain_id(line[1:]) + "\n")
                n += 1
    return f"chain list: {n} polypeptide(L) chains -> {out_path}"


def tsv_key_list(tsv_path: str | Path, out_path: str | Path) -> str:
    """Write column-0 keys of a TSV, one per line (e.g. seqid_to_chains -> seq_ids)."""
    n = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(tsv_path).open() as f, Path(out_path).open("w") as out:
        for line in f:
            key = line.split("\t", 1)[0].strip()
            if key:
                out.write(key + "\n")
                n += 1
    return f"key list: {n} keys -> {out_path}"
