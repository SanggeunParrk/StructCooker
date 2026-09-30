"""Resolve shared input locations for downloads and build configuration."""

from __future__ import annotations

import os
from pathlib import Path


def mmcif_root(data_root: str | Path) -> Path:
    """Return the mmCIF input directory, matching db/pdb/cif.yaml."""
    return Path(os.environ.get("MMCIF_ROOT", str(Path(data_root) / "BioMol/materials/raw/cif")))


def distillation_root(data_root: str | Path) -> Path:
    """Return the supplied OpenFold input directory, matching db/distillation."""
    return Path(os.environ.get("DISTILLATION_ROOT", str(Path(data_root) / "BioMol/materials/raw/openfold_distillation")))
