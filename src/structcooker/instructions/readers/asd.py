from __future__ import annotations

import hashlib
import re
from pathlib import Path

# chain_set -> (role in ASD's `chains` struct, header tag, Auth letter), in the order the sequences are emitted
_CHAIN_SETS = {
    "antibody": (("antibody_chains", "heavy", "H"), ("antibody_chains", "light", "L")),
    "antigen": (("target_chains", "full", "T"),),
}
_AUTH = {"H": "H", "L": "L", "T": "A"}
_LETTERS = re.compile(r"[A-Z]+")


def asd_header(tag: str, sequence: str) -> str:
    """The FASTA header of an ASD chain: ``ASD<H|L|T>_<sha1[:12]> | polypeptide(L) | Auth:<H|L|A>``."""
    digest = hashlib.sha1(sequence.encode()).hexdigest()[:12]  # noqa: S324 (a name, not a secret)
    return f"ASD{tag}_{digest} | polypeptide(L) | Auth:{_AUTH[tag]}"


def load_asd_fasta(asd_dir: str | Path, chain_set: str) -> dict[str, str]:
    """One FASTA of the distinct chain sequences of the ASD Delta table, as ``{header: sequence}``.

    ``chain_set`` is ``antibody`` (heavy and light chains; a VHH or scFv is a ``heavy``) or
    ``antigen`` (``chains.target_chains.full``). Sequences are upper-cased and stripped; an empty chain is
    skipped, and so is a sequence with a character outside A-Z (v2.0 has six, all in ``literature``: ``*``, ``;``,
    ``(``, ``)``, ``?`` inside the chain; counted and shown in the log). Residues outside the 20 standard ones
    (X, B, Z, U, O) are kept as they are: the id names the sequence ASD ships.

    The header is ``ASD<tag>_<sha1[:12]> | polypeptide(L) | Auth:<H|L|A>`` -- the layout the seq_id map and
    the clustering read. The order is fixed (role, length, sequence) because a seq_id is a counter assigned
    in order of first appearance: the same table must give the same ids.
    """
    import pyarrow.parquet as pq  # not needed by any other reader

    asd_dir = Path(asd_dir)
    files = sorted(asd_dir.glob("part-*.parquet"))
    if not files:
        msg = f"no part-*.parquet under {asd_dir}"
        raise FileNotFoundError(msg)
    per_role: dict[str, set[str]] = {tag: set() for _, _, tag in _CHAIN_SETS[chain_set]}
    n_rows = 0
    n_empty = dict.fromkeys(per_role, 0)
    invalid: list[str] = []
    for path in files:
        for row in pq.read_table(path, columns=["chains"]).column("chains").to_pylist():
            n_rows += 1
            for field, key, tag in _CHAIN_SETS[chain_set]:
                raw = dict(row[field] or []).get(key)
                sequence = "".join((raw or "").split()).upper()
                if not sequence:
                    n_empty[tag] += 1
                    continue
                if not _LETTERS.fullmatch(sequence):
                    invalid.append(f"{key}:{sequence[:50]}")
                    continue
                per_role[tag].add(sequence)

    fasta: dict[str, str] = {}
    for _, _, tag in _CHAIN_SETS[chain_set]:
        for sequence in sorted(per_role[tag], key=lambda s: (len(s), s)):
            header = asd_header(tag, sequence)
            if fasta.setdefault(header, sequence) != sequence:
                msg = f"header collision: {header}"
                raise ValueError(msg)
    print(  # noqa: T201 (CLI progress)
        f"asd {chain_set}: {n_rows} rows, "
        + ", ".join(
            f"{tag}: {len(per_role[tag])} distinct, {n_empty[tag]} rows without"
            for tag in per_role
        )
        + f"; {len(invalid)} rows skipped for a character outside A-Z {invalid[:6]}",
    )
    return fasta
