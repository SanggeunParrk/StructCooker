"""The one root everything lives under, and the shared external search databases.

``BIOMOL_ROOT`` (default ``/data/shared/cssb_data/BioMol_clean``) holds the whole release:
``materials/raw`` (downloads as fetched), ``materials/reference`` (supplied inputs),
``materials/intermediate`` (derived build inputs), ``lmdb`` (the DBs), ``metadata`` (released
tables) and ``docs``. Every db config resolves its paths against it, so moving or renaming
the tree is one variable. ``SEQ_DB_ROOT`` points at the shared sequence search databases
(UniRef30, MMseqs2), which are not part of the release.
"""
from __future__ import annotations

import os
from pathlib import Path

DEFAULT_ROOT = "/data/shared/cssb_data/BioMol_clean"
DEFAULT_SEQ_DB_ROOT = "/data/shared/cssb_data/seq_db"


def biomol_root() -> Path:
    """The release tree (``BIOMOL_ROOT``)."""
    return Path(os.environ.get("BIOMOL_ROOT", DEFAULT_ROOT))


def seq_db_root() -> Path:
    """The shared sequence search databases (``SEQ_DB_ROOT``)."""
    return Path(os.environ.get("SEQ_DB_ROOT", DEFAULT_SEQ_DB_ROOT))


def mmcif_root(root: str | Path) -> Path:
    """The mmCIF mirror, matching db/pdb/cif.yaml."""
    return Path(root) / "materials/raw/cif"


def distillation_root(root: str | Path) -> Path:
    """The OpenFold distillation download, matching db/distillation."""
    return Path(root) / "materials/raw/openfold_distillation"
