"""Antibody / TCR chain classification from sequence (ANARCI, IMGT numbering).

Returns the ANARCI chain type: ``H`` (heavy), ``K`` / ``L`` (light), ``A`` / ``B`` / ``G`` / ``D``
(TCR alpha / beta / gamma / delta), or ``None``. Nanobodies and scFvs come out as one chain.
Results are cached per sequence, optionally in a TSV (``seq<TAB>type``) so a full build can
classify every unique sequence once (:func:`annotate_antibodies`, the ``unit/antibody_cache``
stage).
"""

from __future__ import annotations

import multiprocessing as mp
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable

# complementary chain types that are always bundled
PAIR_PARTNERS = {
    "H": ("K", "L"),
    "K": ("H",),
    "L": ("H",),
    "A": ("B",),
    "B": ("A",),
    "G": ("D",),
    "D": ("G",),
}


class AntibodyClassifier:
    """Classify sequences as antibody / TCR chains with ANARCI, with a per-sequence cache."""

    def __init__(
        self,
        cache_tsv: str | Path | None = None,
        min_len: int = 70,
        max_len: int = 600,
        lazy: bool = True,
    ) -> None:
        self.cache: dict[str, str | None] = {}
        self.min_len, self.max_len, self.lazy = min_len, max_len, lazy
        if cache_tsv and Path(cache_tsv).exists():
            for line in Path(cache_tsv).read_text().splitlines():
                s, _, t = line.partition("\t")
                self.cache[s] = t or None

    def __call__(self, seq: str) -> str | None:
        """Return the ANARCI chain type of ``seq`` (cached), or None."""
        if seq in self.cache:
            return self.cache[seq]
        res = None
        if not self.lazy:  # frozen cache: an unknown sequence is not an antibody
            return None
        if self.min_len <= len(seq) <= self.max_len:
            from anarci import anarci  # imported lazily: hmmer must be on PATH

            num, det, _ = anarci([("q", seq)], scheme="imgt", output=False)
            if num[0] and det[0]:
                res = det[0][0].get("chain_type")
        self.cache[seq] = res
        return res

    def save(self, path: str | Path) -> None:
        """Write the cache as a ``seq<TAB>type`` TSV."""
        Path(path).write_text(
            "".join(f"{s}\t{t or ''}\n" for s, t in self.cache.items())
        )


def _classify_chunk(chunk: list[str]) -> list[tuple[str, str | None]]:
    """ANARCI chain type of every sequence of ``chunk``; one bad sequence must not sink the chunk."""
    from anarci import anarci

    def one(batch: list[str]) -> list[str | None]:
        num, det, _ = anarci(
            [(f"s{i}", s) for i, s in enumerate(batch)], scheme="imgt", output=False
        )
        return [
            (d[0].get("chain_type") if (n and d) else None)
            for n, d in zip(num, det, strict=True)
        ]

    def alone(seq: str) -> str | None:
        try:
            return one([seq])[0]
        except Exception:  # noqa: BLE001  (an unclassifiable sequence is "not an antibody")
            return None

    try:
        types = one(chunk)
    except Exception:  # noqa: BLE001  (a failing batch is retried sequence by sequence)
        types = [alone(s) for s in chunk]
    return list(zip(chunk, types, strict=True))


def unique_sequences(fasta_path: str | Path, min_len: int, max_len: int) -> list[str]:
    """Sorted unique sequences of a one-line-per-sequence FASTA with ``min_len <= length <= max_len``."""
    seqs: set[str] = set()
    with Path(fasta_path).open() as handle:
        for line in handle:
            if not line.startswith(">"):
                s = line.strip()
                if min_len <= len(s) <= max_len:
                    seqs.add(s)
    return sorted(seqs)


def annotate_antibodies(
    fasta_path: str | Path,
    out_path: str | Path,
    n_jobs: int = 1,
    chunk_size: int = 200,
    classify: Callable[[list[str]], list[tuple[str, str | None]]] = _classify_chunk,
) -> str:
    """Classify every unique sequence of ``fasta_path`` with ANARCI and write ``out_path``.

    Output lines are ``sequence<TAB>type`` in sequence order (an empty type is "not an antibody
    or TCR"), the cache :class:`AntibodyClassifier` reads. Only lengths 70-600 are sent to ANARCI.
    Returns a one-line status.
    """
    seqs = unique_sequences(fasta_path, 70, 600)
    chunks = [seqs[i : i + chunk_size] for i in range(0, len(seqs), chunk_size)]
    types: dict[str, str | None] = {}
    if n_jobs > 1:
        with mp.get_context("fork").Pool(n_jobs) as pool:
            for result in pool.imap_unordered(classify, chunks):
                types.update(result)
    else:
        for chunk in chunks:
            types.update(classify(chunk))
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".part")
    tmp.write_text("".join(f"{s}\t{types[s] or ''}\n" for s in seqs))
    tmp.replace(out)
    n_ab = sum(1 for s in seqs if types[s])
    return f"{n_ab} of {len(seqs)} sequences classified as antibody/TCR"
