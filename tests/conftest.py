"""Shared pytest fixtures for StructCooker tests."""

from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path

import pytest

FIXTURE_DIR = Path(__file__).parent / "fixtures"
DB_AVAILABLE = find_spec("lmdb") is not None and find_spec("joblib") is not None


@pytest.fixture(scope="session")
def fixture_dir() -> Path:
    """Directory holding small per-component CCD `.cif` fixtures."""
    return FIXTURE_DIR


@pytest.fixture
def built_ccd_lmdb(tmp_path: Path) -> Path:
    """Build a tiny CCD LMDB from the bundled fixtures and return its path.

    Uses the real ingest recipe + reader/writer (the same ones
    `configs/ingest/ccd_lmdb.yaml` wires), so this exercises build -> store
    end-to-end on a handful of components.
    """
    if not DB_AVAILABLE:
        pytest.skip("lmdb/joblib not installed")

    from datacooker.lmdb import build_lmdb

    import structcooker.workflows.ingest.ccd as ccd_recipe_module
    from structcooker.instructions.readers.cif import get_cif_data
    from structcooker.instructions.transforms.codecs import to_bytes

    cifs = sorted(FIXTURE_DIR.glob("*.cif"))
    env_path = tmp_path / "ccd_fixtures.lmdb"
    build_lmdb(
        *cifs,
        env_path=env_path,
        recipe=Path(ccd_recipe_module.__file__),
        loader=get_cif_data,
        serializer=to_bytes,
        n_jobs=1,
        test_run=False,
    )
    return env_path
