"""Antibody / TCR chain classification from sequence (ANARCI, IMGT numbering).

Returns the ANARCI chain type: ``H`` (heavy), ``K`` / ``L`` (light), ``A`` / ``B`` / ``G`` / ``D``
(TCR alpha / beta / gamma / delta), or ``None``. Nanobodies and scFvs come out as one chain.
Results are cached per sequence, optionally in a TSV (``seq<TAB>type``) so a full build can
classify every unique sequence once (see ``scripts/interface/annotate_antibodies.py``).
"""

from __future__ import annotations

from pathlib import Path

# complementary chain types that are always bundled
PAIR_PARTNERS = {"H": ("K", "L"), "K": ("H",), "L": ("H",), "A": ("B",), "B": ("A",), "G": ("D",), "D": ("G",)}


class AntibodyClassifier:
    def __init__(self, cache_tsv: str | Path | None = None, min_len: int = 70, max_len: int = 600, lazy: bool = True):
        self.cache: dict[str, str | None] = {}
        self.min_len, self.max_len, self.lazy = min_len, max_len, lazy
        if cache_tsv and Path(cache_tsv).exists():
            for line in Path(cache_tsv).read_text().splitlines():
                s, _, t = line.partition("\t")
                self.cache[s] = t or None

    def __call__(self, seq: str) -> str | None:
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
        Path(path).write_text("".join(f"{s}\t{t or ''}\n" for s, t in self.cache.items()))
