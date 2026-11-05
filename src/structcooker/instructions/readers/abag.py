from __future__ import annotations

import csv
import hashlib
import re
from collections.abc import Iterator
from pathlib import Path
from typing import NamedTuple

from .abag_binding import abrank_verdict, poisk_verdict
from .asd_binding import binding_verdict

_LETTERS = re.compile(r"[A-Z]+")
_AUTH = {"H": "H", "L": "L", "T": "A"}
# the FR4 end of a V domain: heavy J segment / light J segments (psh's sanity check on cssb-master, docs/PAIRS_EXPORT.md 7-1)
_HEAVY_FR4 = re.compile(r"WG.G[TA]|VTVSS")
_LIGHT_FR4 = re.compile(r"FG.GT|[KV][LVTE][ED]I[KR]$|VTVL$")


class Row(NamedTuple):
    """One measurement row of a source: the chains, what it says about binding (``True`` / ``False`` / ``None`` no verdict) and where it comes from."""

    source: str
    dataset: str
    heavy: str
    light: str
    antigen: str
    verdict: bool | None
    nanobody: bool
    scfv: bool
    name: str


def normalize(sequence: str | None) -> str:
    return "".join((sequence or "").split()).upper()


def abag_header(tag: str, sequence: str) -> str:
    """The FASTA header of a chain of the antibody-antigen set: ``ABAG<H|L|T>_<sha1[:12]> | polypeptide(L) | Auth:<H|L|A>``."""
    digest = hashlib.sha1(sequence.encode()).hexdigest()[:12]  # noqa: S324 (a name, not a secret)
    return f"ABAG{tag}_{digest} | polypeptide(L) | Auth:{_AUTH[tag]}"


def columns_swapped(heavy: str, light: str) -> bool:
    """AbRank's heavy and light columns are swapped in 807 rows: the heavy column ends in a light FR4 *and* the light column in a heavy FR4 (both, so a
    rare J variant of one chain is not taken for a swap).
    """
    return bool(_LIGHT_FR4.search(heavy[-20:]) and _HEAVY_FR4.search(light[-20:]))


def iter_asd_rows(asd_dir: str | Path) -> Iterator[Row]:
    import pyarrow.parquet as pq  # not needed by any other reader

    files = sorted(Path(asd_dir).glob("part-*.parquet"))
    if not files:
        msg = f"no part-*.parquet under {asd_dir}"
        raise FileNotFoundError(msg)
    for path in files:
        columns = [
            "dataset",
            "affinity",
            "affinity_type",
            "chains",
            "metadata",
            "nanobody",
            "scfv",
        ]
        for row in pq.read_table(path, columns=columns).to_pylist():
            antibody = dict(row["chains"]["antibody_chains"] or [])
            target = dict(row["chains"]["target_chains"] or [])
            yield Row(
                "asd",
                row["dataset"],
                normalize(antibody.get("heavy")),
                normalize(antibody.get("light")),
                normalize(target.get("full")),
                binding_verdict(row["affinity_type"], row["affinity"], row["dataset"]),
                bool(row["nanobody"]),
                bool(row["scfv"]),
                (row["metadata"] or {}).get("target_name") or "",
            )


def iter_abrank_rows(path: str | Path) -> Iterator[Row]:
    """AbRank (a tab-separated table). An antigen of several chains (``:`` between them) is not a single chain and is left out."""
    csv.field_size_limit(10**9)
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            antigen = normalize(row["Ag_seq"])
            if ":" in antigen:
                continue
            heavy, light = (
                normalize(row["Ab_heavy_chain_seq"]),
                normalize(row["Ab_light_chain_seq"]),
            )
            if columns_swapped(heavy, light):
                heavy, light = light, heavy
            verdict = abrank_verdict(
                row["Affinity_Kd [nM]"], row["IC50 [ug/mL]"], row["Aff_op"] or "="
            )
            yield Row(
                "abrank",
                row["Source"],
                heavy,
                light,
                antigen,
                verdict,
                False,
                False,
                row["Ag_name"],
            )


def iter_poisk_rows(path: str | Path) -> Iterator[Row]:
    csv.field_size_limit(10**9)
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle):
            keep, verdict = poisk_verdict(
                row["target.binding_affinity.affinity_value"], row["audit.verdict"]
            )
            if not keep:
                continue
            yield Row(
                "poisk",
                row["antibody_type"],
                normalize(row["vh.sequence_from_patent"]),
                normalize(row["vl.sequence_from_patent"]),
                normalize(row["target.sequence_from_patent"]),
                verdict,
                row["antibody_type"] == "VHH",
                False,
                row["target.protein_name"],
            )


def iter_rows(
    asd_dir: str | Path, abrank_path: str | Path, poisk_path: str | Path
) -> Iterator[Row]:
    yield from iter_asd_rows(asd_dir)
    yield from iter_abrank_rows(abrank_path)
    yield from iter_poisk_rows(poisk_path)


def valid(row: Row) -> bool:
    """A usable row: a heavy chain and an antigen, every chain A-Z only."""
    return bool(
        row.heavy
        and row.antigen
        and _LETTERS.fullmatch(row.heavy)
        and _LETTERS.fullmatch(row.antigen)
        and (not row.light or _LETTERS.fullmatch(row.light)),
    )


def load_abag_fasta(
    asd_dir: str | Path, abrank_path: str | Path, poisk_path: str | Path, chain_set: str
) -> dict[str, str]:
    """One FASTA of the distinct chain sequences of ASD, AbRank and POISK (every usable row): ``chain_set`` is ``antigen`` or ``antibody`` (heavy, then light).

    A sequence several sources hold is one record. The order is (role, length, sequence): a seq_id is a counter assigned in order of first appearance, so
    the same inputs must give the same ids.
    """
    per_role: dict[str, set[str]] = {"T": set(), "H": set(), "L": set()}
    for row in iter_rows(asd_dir, abrank_path, poisk_path):
        if not valid(row):
            continue
        per_role["T"].add(row.antigen)
        per_role["H"].add(row.heavy)
        if row.light:
            per_role["L"].add(row.light)
    fasta: dict[str, str] = {}
    for tag in ("T",) if chain_set == "antigen" else ("H", "L"):
        for sequence in sorted(per_role[tag], key=lambda s: (len(s), s)):
            fasta.setdefault(abag_header(tag, sequence), sequence)
    return fasta
