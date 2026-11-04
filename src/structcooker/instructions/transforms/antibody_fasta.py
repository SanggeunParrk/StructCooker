"""The antibody chains of a PDB fasta, and the CDR regions of those chains.

Two fasta-to-fasta steps (recipes ``exports/ab_seq.py`` and ``exports/cdr_region.py``):

* :func:`filter_fasta_antibody` keeps the records of an antibody sequence, in the format of the input
  (``>100D_A_. | polypeptide(L) | Auth:A``). "Antibody" is the rule ``separate_sequences`` clusters with: a
  sequence is an antibody if it is the heavy or light chain of a SAbDab entry (matched on the PDB id and the
  author chain id), and then every record carrying that sequence is one, SAbDab-listed or not.
* :func:`extract_cdr_regions` numbers each distinct antibody sequence with ANARCI (Chothia) and writes its CDR
  residues, one record per chain and CDR: the header of the chain plus ``| CDR:H3``.

The CDR extents are ``sequence.chotia_map`` (L1 24-34, L2 50-56, L3 89-97, H1 26-32, H2 50-65, H3 95-102): insertion
codes count with their number, so H3 includes 100A, 100B... The H2 extent is the Kabat one, wider than Chothia's
52-56. ``ranges`` overrides them.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

from anarci import run_anarci

from structcooker.instructions.transforms.sequence import chotia_map

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping, Sequence

_L_TYPE = "polypeptide(L)"
CDR_NAMES = ("H1", "H2", "H3", "L1", "L2", "L3")
# ANARCI chain types: H is a heavy chain, K and L are the two light-chain classes.
_REGION_OF_CHAIN_TYPE = {"H": "H", "K": "L", "L": "L"}


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """``(header, sequence)`` of every record; ``header`` is the line without ``>``."""
    header, parts = None, []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].rstrip("\n"), []
            elif header is not None:
                parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)


def _is_polypeptide_l(header: str) -> bool:
    fields = header.split("|")
    return len(fields) > 1 and fields[1].strip() == _L_TYPE


def _auth_key(header: str) -> str | None:
    """``<PDB id>_<author chain id>`` of a record, the key SAbDab's chain columns are matched on."""
    if "Auth:" not in header:
        return None
    return header.split("|")[0].strip().split("_")[0] + "_" + header.split("Auth:")[1].strip()


def sabdab_antibody_sequences(fasta_path: str | Path, sabdab_summary_path: str | Path) -> tuple[set[str], set[str]]:
    """Heavy-chain and light-chain sequences of the fasta that SAbDab lists (the ``separate_sequences`` rule)."""
    heavy_keys: set[str] = set()
    light_keys: set[str] = set()
    with Path(sabdab_summary_path).open(encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            items = line.rstrip("\n").split("\t")
            pdb_id, h_chain, l_chain = items[0], items[1], items[2]
            if h_chain != "NA":
                heavy_keys.add(f"{pdb_id}_{h_chain}".upper())
            if l_chain != "NA":
                light_keys.add(f"{pdb_id}_{l_chain}".upper())
    wanted = heavy_keys | light_keys
    seq_of: dict[str, str] = {}
    for header, sequence in read_fasta(fasta_path):
        if _is_polypeptide_l(header):
            key = _auth_key(header)
            if key in wanted:
                seq_of[key] = sequence          # as separate_sequences: the last record of an author chain
    heavy = {seq_of[k] for k in heavy_keys if k in seq_of}
    light = {seq_of[k] for k in light_keys if k in seq_of}
    return heavy, light


def filter_fasta_antibody(fasta_path: str | Path, sabdab_summary_path: str | Path, out_path: str | Path) -> str:
    """Write the records of ``fasta_path`` whose sequence is an antibody heavy or light chain, headers unchanged."""
    heavy, light = sabdab_antibody_sequences(fasta_path, sabdab_summary_path)
    antibody = heavy | light
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with Path(out_path).open("w", encoding="utf-8") as out:
        for header, sequence in read_fasta(fasta_path):
            if sequence in antibody and _is_polypeptide_l(header):
                out.write(f">{header}\n{sequence}\n")
                kept += 1
    return (f"antibody fasta: {kept} records, {len(antibody)} distinct sequences "
            f"({len(heavy)} heavy, {len(light)} light, {len(heavy & light)} both) -> {out_path}")


def _numbered_domains(
    sequences: Sequence[str],
    ncpu: int,
    scheme: str,
) -> list[list[tuple[str, list[tuple[int, str]]]]]:
    """Per sequence, its antibody domains as ``(H | L, [(Chothia number, residue), ...])``; ``[]`` when ANARCI finds none."""
    queries = [(str(i), s.replace("X", "")) for i, s in enumerate(sequences) if s.replace("X", "")]
    if not queries:
        return [[] for _ in sequences]
    _, numbering, details = run_anarci(queries, scheme=scheme, ncpu=ncpu, allow={"H", "K", "L"})[:3]
    found: dict[int, list[tuple[str, list[tuple[int, str]]]]] = {}
    for (name, _), domains, info in zip(queries, numbering, details, strict=True):
        if not domains:
            continue
        found[int(name)] = [
            (_REGION_OF_CHAIN_TYPE[meta["chain_type"]], [(idx, aa) for (idx, _ins), aa in domain[0]])
            for domain, meta in zip(domains, info, strict=True)
        ]
    return [found.get(i, []) for i in range(len(sequences))]


def cdr_regions_of(
    domains: Sequence[tuple[str, Sequence[tuple[int, str]]]],
    cdrs: Sequence[str] = CDR_NAMES,
    ranges: Mapping[str, Sequence[int]] | None = None,
) -> list[tuple[str, str]]:
    """``[(cdr name, residues)]`` of numbered domains; a second domain of the same type is named ``H3.2``.

    A CDR with no residue in the numbered sequence (a truncated chain) is left out.
    """
    ranges = chotia_map if ranges is None else ranges
    seen: dict[str, int] = defaultdict(int)
    regions = []
    for region, residues in domains:
        seen[region] += 1
        suffix = "" if seen[region] == 1 else f".{seen[region]}"
        for cdr in cdrs:
            if not cdr.startswith(region):
                continue
            extent = set(ranges[cdr])
            residue_string = "".join(aa for idx, aa in residues if idx in extent and aa != "-")
            if residue_string:
                regions.append((cdr + suffix, residue_string))
    return regions


def extract_cdr_regions(
    fasta_path: str | Path,
    out_path: str | Path,
    cdrs: Sequence[str] = CDR_NAMES,
    scheme: str = "chothia",
    batch_size: int = 2000,
    ncpu: int = 8,
    ranges: Mapping[str, Sequence[int]] | None = None,
    failed_path: str | Path | None = None,
) -> str:
    """CDR residues of every antibody record, one fasta record per chain and CDR.

    Each distinct sequence is numbered once. Records with no antibody domain (ANARCI finds none) are listed in
    ``failed_path``, one header per line, and get no CDR record.
    """
    unknown = [c for c in cdrs if c not in CDR_NAMES]
    if unknown:
        msg = f"unknown CDR names {unknown}; known: {list(CDR_NAMES)}"
        raise ValueError(msg)
    records = list(read_fasta(fasta_path))
    distinct = sorted({sequence for _, sequence in records})
    regions_of: dict[str, list[tuple[str, str]]] = {}
    for start in range(0, len(distinct), batch_size):
        chunk = distinct[start : start + batch_size]
        for sequence, domains in zip(chunk, _numbered_domains(chunk, ncpu, scheme), strict=True):
            regions_of[sequence] = cdr_regions_of(domains, cdrs, ranges)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    n_out = 0
    failed = []
    with Path(out_path).open("w", encoding="utf-8") as out:
        for header, sequence in records:
            regions = regions_of[sequence]
            if not regions:
                failed.append(header)
            for name, residues in regions:
                out.write(f">{header} | CDR:{name}\n{residues}\n")
                n_out += 1
    if failed_path is not None:
        Path(failed_path).write_text("".join(f"{h}\n" for h in failed), encoding="utf-8")
    return (f"CDR fasta: {n_out} records from {len(records) - len(failed)} of {len(records)} chains "
            f"({len(distinct)} distinct sequences, {len(failed)} chains without a CDR) -> {out_path}")
