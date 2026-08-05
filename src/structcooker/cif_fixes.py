"""Apply the manual mmCIF substitutions the pdb/cif build needs.

A small set of PDB entries error out (or build wrongly) from the *current* wwPDB
mmCIF snapshot -- mostly NMR ensembles whose non-polymer ligand/ion is re-numbered
per model, which breaks the atom->scheme match. The original BioMol build worked
around this by replacing each entry's mmCIF with a version from an older, known-good
snapshot before ingest (legacy ``scripts/manually_fix_cif.py``). This module ports that
step: it overlays the corrected cif for each listed id into the mmCIF input directory.

The id list lives in ``db/pdb/manual_cif_fixes.txt``; the corrected cifs are a *provided*
external input (the old snapshot is not redistributed here). See docs/manual-cif-fixes.md.
"""
from __future__ import annotations

import gzip
import shutil
from pathlib import Path

# Written into the mmCIF dir after a successful apply, so ``inspect`` can tell whether
# the substitutions are in place without re-reading every corrected file.
APPLIED_MARKER = ".manual_cif_fixes_applied"


def load_fix_ids(list_path: Path) -> list[str]:
    """Return the PDB ids from ``manual_cif_fixes.txt`` (skips blanks and ``#`` lines)."""
    ids: list[str] = []
    for raw in list_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            ids.append(line.lower())
    return ids


def _source_cif(source_dir: Path, pdb_id: str) -> Path | None:
    """Locate ``pdb_id``'s corrected cif under ``source_dir`` (divided or flat, gz or not)."""
    sub = pdb_id[1:3]
    for candidate in (
        source_dir / sub / f"{pdb_id}.cif.gz",
        source_dir / sub / f"{pdb_id}.cif",
        source_dir / f"{pdb_id}.cif.gz",
        source_dir / f"{pdb_id}.cif",
    ):
        if candidate.exists():
            return candidate
    return None


def apply_fixes(
    ids: list[str],
    source_dir: Path,
    mmcif_dir: Path,
    *,
    dry_run: bool = False,
) -> tuple[list[str], list[str]]:
    """Overlay each id's corrected cif from ``source_dir`` into ``mmcif_dir``.

    Returns ``(applied, missing)`` -- ids written, and ids with no corrected cif in the
    source. Gzipped sources are decompressed to ``<id>.cif``. Writes an
    :data:`APPLIED_MARKER` listing the applied ids so a later ``inspect`` can see them.
    """
    applied: list[str] = []
    missing: list[str] = []
    for pdb_id in ids:
        src = _source_cif(source_dir, pdb_id)
        if src is None:
            missing.append(pdb_id)
            continue
        if dry_run:
            applied.append(pdb_id)
            continue
        dst = mmcif_dir / f"{pdb_id}.cif"
        if src.suffix == ".gz":
            with gzip.open(src, "rb") as fh_in, dst.open("wb") as fh_out:
                shutil.copyfileobj(fh_in, fh_out)
        else:
            shutil.copyfile(src, dst)
        applied.append(pdb_id)
    if applied and not dry_run:
        (mmcif_dir / APPLIED_MARKER).write_text(
            "\n".join(sorted(applied)) + "\n", encoding="utf-8",
        )
    return applied, missing


def applied_ids(mmcif_dir: Path) -> set[str]:
    """Return the ids recorded as applied in ``mmcif_dir`` (empty if none / no marker)."""
    marker = mmcif_dir / APPLIED_MARKER
    if not marker.exists():
        return set()
    return {
        line.strip().lower()
        for line in marker.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
