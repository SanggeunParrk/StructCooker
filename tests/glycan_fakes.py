"""Synthetic BioMol-shaped values for the Re-Glyco tests (shared by tests/glycan)."""

from __future__ import annotations

import numpy as np

from structcooker.instructions.transforms.delta_codec import RecordRef

PDB = "2pmv"


def value_for(key: str, n_protein: int = 4000) -> dict:
    """Build a BioMol-shaped value: one protein context for every record, a glycan that depends on the record."""
    rng = np.random.default_rng(
        0
    )  # the protein is the same in every record of the entry
    protein = rng.normal(size=(n_protein, 3)).astype(np.float32)
    names = np.array([f"C{i % 17}" for i in range(n_protein)])
    site, glycan, conf = key.split("_")[1], key.split("_")[2], int(key.split("conf")[1])
    n_glycan = 20 + 3 * int(
        glycan.removeprefix("GS")[-1]
    )  # another glycan has another atom count
    g = (
        np.random.default_rng(1000 * int(glycan.removeprefix("GS")) + conf)
        .normal(size=(n_glycan, 3))
        .astype(np.float32)
    )
    return {
        "assembly_dict": {
            "1_1_.": {
                "atoms": {
                    "xyz": np.concatenate([protein, g]),
                    "name": np.concatenate([names, np.array(["O"] * n_glycan)]),
                },
                "index": np.arange(n_protein + n_glycan, dtype=np.int32),
            }
        },
        "metadata_dict": {
            "id": key,
            "deposition_date": "2020-01-01",
            "site": site,
            "resolution": 2.1,
        },
    }


def make_refs(spec: dict[str, int]) -> list[RecordRef]:
    """``{group: n conformers}`` -> the refs of one entry; group is ``<site>_<glycan>``."""
    return [
        RecordRef(f"{PDB}_{group}_conf{c:03d}", group, c)
        for group, n in spec.items()
        for c in range(n)
    ]


SPEC = {"A293_GS00001": 3, "A293_GS00007": 4, "B12_GS00003": 2}


def same(a, b):
    """Equality of two decoded values, arrays by dtype and content."""
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, np.ndarray):
        return a.dtype == b.dtype and np.array_equal(a, b)
    return a == b
