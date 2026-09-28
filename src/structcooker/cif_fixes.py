"""Apply the supplied substitutions for a historical BioMol mmCIF snapshot.

The optional ID list lives in ``db/pdb/manual_cif_fixes.txt``. This is not a
universal preprocessing requirement for new releases; see docs/manual-cif-fixes.md.
"""
from __future__ import annotations

import gzip
import shutil
import tempfile
from pathlib import Path

# Bookkeeping for historical substitutions; this is not a content checksum.
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
    source. Replaces the existing divided/flat ``<id>.cif.gz`` read by pdb/cif.
    Each replacement is staged before an atomic rename. Writes an
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
        candidates = [p for p in (
            mmcif_dir / pdb_id[1:3] / f"{pdb_id}.cif.gz",
            mmcif_dir / f"{pdb_id}.cif.gz",
        ) if p.exists()]
        if len(candidates) > 1:
            msg = f"Duplicate input files for {pdb_id}: {candidates}"
            raise ValueError(msg)
        dst = candidates[0] if candidates else mmcif_dir / f"{pdb_id}.cif.gz"
        dst.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=dst.parent, delete=False) as staged:
            stage_path = Path(staged.name)
        try:
            opener = gzip.open if src.suffix == ".gz" else Path.open
            with opener(src, "rb") as fh_in, gzip.open(stage_path, "wb") as fh_out:
                shutil.copyfileobj(fh_in, fh_out)
            stage_path.replace(dst)
        finally:
            stage_path.unlink(missing_ok=True)
        applied.append(pdb_id)
    if applied and not dry_run:
        (mmcif_dir / APPLIED_MARKER).write_text(
            "\n".join(sorted(applied_ids(mmcif_dir) | set(applied))) + "\n", encoding="utf-8",
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
