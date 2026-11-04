"""The antibody chains of a PDB fasta, and the CDR positions of their sequences.

Two steps (recipes ``exports/ab_seq.py`` and ``exports/cdr_region.py``):

* :func:`filter_fasta_antibody` keeps the records of an antibody sequence, in the format of the input
  (``>100D_A_. | polypeptide(L) | Auth:A``). "Antibody" is the rule ``separate_sequences`` clusters with: a
  sequence is an antibody if it is the heavy or light chain of a SAbDab entry (matched on the PDB id and the
  author chain id), and then every record carrying that sequence is one, SAbDab-listed or not.
* :func:`extract_cdr_table` numbers each distinct antibody sequence with ANARCI (Chothia) and writes a CSV, one row
  per sequence: ``seqid`` and the start and end of CDR1, CDR2 and CDR3. CDRn is H1/H2/H3 of a heavy domain and
  L1/L2/L3 of a light one.

Positions are 0-based, the end exclusive, on the sequence exactly as the fasta and ``seq_id_map`` hold it
(``X`` residues included): the CDR is ``sequence[start:end]``.

The CDR extents are ``sequence.chotia_map`` (L1 24-34, L2 50-56, L3 89-97, H1 26-32, H2 50-65, H3 95-102): insertion
codes count with their number, so H3 includes 100A, 100B... The H2 extent is the Kabat one, wider than Chothia's
52-56. ``ranges`` overrides them.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import TYPE_CHECKING

from anarci import run_anarci

from structcooker.instructions.transforms.sequence import chotia_map

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping, Sequence

_L_TYPE = "polypeptide(L)"
CDR_COLUMNS = [f"cdr{n}_{edge}" for n in (1, 2, 3) for edge in ("start", "end")]
# ANARCI chain types: H is a heavy chain, K and L are the two light-chain classes.
_REGION_OF_CHAIN_TYPE = {"H": "H", "K": "L", "L": "L"}

# One numbered domain: its region (H or L) and, for each numbered residue, (Chothia number, position in the sequence).
Domain = tuple[str, list[tuple[int, int]]]


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


def _numbered_domains(sequences: Sequence[str], ncpu: int, scheme: str) -> list[list[Domain]]:
    """Per sequence, its antibody domains with each numbered residue's position; ``[]`` when ANARCI finds none.

    ``X`` residues are not sent to ANARCI; positions are mapped back, so they index the sequence as given.
    """
    kept_of = [[j for j, c in enumerate(s) if c != "X"] for s in sequences]
    queries = [(str(i), "".join(sequences[i][j] for j in kept)) for i, kept in enumerate(kept_of) if kept]
    if not queries:
        return [[] for _ in sequences]
    _, numbering, details = run_anarci(queries, scheme=scheme, ncpu=ncpu, allow={"H", "K", "L"})[:3]
    found: dict[int, list[Domain]] = {}
    for (name, query), domains, info in zip(queries, numbering, details, strict=True):
        if not domains:
            continue
        i, kept, parsed = int(name), kept_of[int(name)], []
        for (domain_numbering, start, _end), meta in zip(domains, info, strict=True):
            at, residues = start, []
            for (idx, _ins), aa in domain_numbering:
                if aa == "-":
                    continue
                if query[at] != aa:
                    msg = f"ANARCI residue {aa} does not sit at position {at} of {query}"
                    raise ValueError(msg)
                residues.append((idx, kept[at]))
                at += 1
            parsed.append((_REGION_OF_CHAIN_TYPE[meta["chain_type"]], residues))
        found[i] = parsed
    return [found.get(i, []) for i in range(len(sequences))]


def cdr_spans_of(domain: Domain, ranges: Mapping[str, Sequence[int]] | None = None) -> list[tuple[int, int] | None]:
    """``(start, end)`` of CDR1, CDR2, CDR3 of a domain (0-based, end exclusive); ``None`` for a CDR with no residue."""
    ranges = chotia_map if ranges is None else ranges
    region, residues = domain
    spans: list[tuple[int, int] | None] = []
    for n in (1, 2, 3):
        extent = set(ranges[f"{region}{n}"])
        positions = [pos for idx, pos in residues if idx in extent]
        spans.append((min(positions), max(positions) + 1) if positions else None)
    return spans


def _row(seq_id: str, spans: Sequence[tuple[int, int] | None]) -> list[str | int]:
    cells: list[str | int] = []
    for span in spans:
        cells += ["", ""] if span is None else list(span)
    return [seq_id, *cells]


def _read_seq_ids(seq_id_map_path: str | Path, sequences: set[str]) -> dict[str, str]:
    """``sequence -> seq_id`` for the polypeptide(L) ids (``P…``) of ``seq_id_map`` that carry one of ``sequences``."""
    ids: dict[str, str] = {}
    with Path(seq_id_map_path).open(encoding="utf-8") as handle:
        for line in handle:
            seq_id, _, sequence = line.rstrip("\n").partition("\t")
            if seq_id.startswith("P") and sequence in sequences:
                ids[sequence] = seq_id
    return ids


def extract_cdr_table(
    fasta_path: str | Path,
    seq_id_map_path: str | Path,
    out_path: str | Path,
    scheme: str = "chothia",
    batch_size: int = 2000,
    ncpu: int = 8,
    ranges: Mapping[str, Sequence[int]] | None = None,
    failed_path: str | Path | None = None,
    extra_path: str | Path | None = None,
) -> str:
    """CSV of the CDR positions of every distinct sequence of ``fasta_path``: ``seqid,cdr1_start,...,cdr3_end``.

    One row per sequence, sorted by seq_id, from its first (N-terminal) antibody domain. The later domains of a
    sequence with several (scFv, tandem) go to ``extra_path`` (``seqid,domain,chain_type`` and the same columns, domain
    2, 3...). Sequences with no antibody domain, or a domain without a residue in any CDR, are listed in ``failed_path`` by
    seq_id and have no row.
    """
    sequences = sorted({sequence for _, sequence in read_fasta(fasta_path)})
    seq_id_of = _read_seq_ids(seq_id_map_path, set(sequences))
    missing = [s for s in sequences if s not in seq_id_of]
    if missing:
        msg = f"{len(missing)} sequences of {fasta_path} have no seq_id in {seq_id_map_path}, e.g. {missing[0][:40]}"
        raise ValueError(msg)
    rows, extra_rows, failed = [], [], []
    for start in range(0, len(sequences), batch_size):
        chunk = sequences[start : start + batch_size]
        for sequence, domains in zip(chunk, _numbered_domains(chunk, ncpu, scheme), strict=True):
            seq_id = seq_id_of[sequence]
            spans = cdr_spans_of(domains[0], ranges) if domains else []
            if all(span is None for span in spans):         # no domain, or one numbered with no residue in any CDR
                failed.append(seq_id)
                continue
            rows.append(_row(seq_id, spans))
            for k, domain in enumerate(domains[1:], start=2):
                extra_rows.append([seq_id, k, domain[0], *_row(seq_id, cdr_spans_of(domain, ranges))[1:]])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(out_path).open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["seqid", *CDR_COLUMNS])
        writer.writerows(sorted(rows))
    if extra_path is not None:
        with Path(extra_path).open("w", encoding="utf-8", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(["seqid", "domain", "chain_type", *CDR_COLUMNS])
            writer.writerows(sorted(extra_rows, key=lambda r: (r[0], r[1])))
    if failed_path is not None:
        Path(failed_path).write_text("".join(f"{s}\n" for s in sorted(failed)), encoding="utf-8")
    return (f"CDR table: {len(rows)} sequences of {len(sequences)} (no usable antibody domain: {len(failed)}; "
            f"{len(extra_rows)} further domains) -> {out_path}")
