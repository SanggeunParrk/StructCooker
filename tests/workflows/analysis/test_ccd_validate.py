"""Integration: validate every item of a freshly built CCD LMDB.

This is the per-item QC the validators exist for — build a small LMDB from the
fixtures with the real ingest recipe, then run the validation *recipe* over each
stored entry and assert the items are internally sound.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from biomol.core.container import FeatureContainer
from biomol.core.feature import EdgeFeature, NodeFeature
from datacooker import execute

from structcooker.workflows.analysis.ccd_validate import RECIPE, TARGETS

pytestmark = pytest.mark.skipif(
    __import__("importlib.util", fromlist=["find_spec"]).find_spec("lmdb") is None,
    reason="lmdb not installed",
)


def _validate(entry: object) -> list[dict]:
    return execute(RECIPE, {"chem_comp": entry}, targets=TARGETS)["issues"]


def test_built_fixtures_are_all_valid(built_ccd_lmdb: Path) -> None:
    from datacooker.lmdb import extract_lmdb_keys, read_lmdb

    from structcooker.instructions.transforms.codecs import from_bytes

    keys = extract_lmdb_keys(built_ccd_lmdb)
    assert set(keys) == {"ala", "hoh", "atp"}

    for key in keys:
        entry = read_lmdb(built_ccd_lmdb, key, deserializer=from_bytes)["chem_comp_dict"]
        issues = _validate(entry)
        assert issues == [], f"{key} unexpectedly flagged: {issues}"


def test_corrupted_entry_is_flagged_through_the_recipe() -> None:
    # A deliberately broken entry (atoms 3.5 A apart but bonded) must be caught
    # when run through the same validation recipe.
    atom = FeatureContainer(features={
        "id": NodeFeature(value=np.array(["C1", "C2"])),
        "element": NodeFeature(value=np.array(["C", "C"])),
        "model_xyz": NodeFeature(value=np.array([["0", "0", "0"], ["3.5", "0", "0"]])),
        "bond_type": EdgeFeature(
            value=np.array(["SING"]),
            src_indices=np.array([0]),
            dst_indices=np.array([1]),
        ),
    })
    issues = _validate({"residue": None, "atom": atom})
    assert any(i["code"] == "bond_too_long" for i in issues)
