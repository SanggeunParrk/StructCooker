"""Fetch the raw external inputs the BioMol pipeline reads.

One function per ``structcooker download`` target. Small public inputs (CCD, SabDab)
download and unpack directly; the huge ones (mmCIF, OpenFold distillation) require an
explicit ``confirmed=True`` and print a size warning first, so nobody kicks off a
multi-hundred-GB / multi-TB transfer by accident. The ``seq_id_map`` seed is NOT here --
it is a Hugging Face dataset (``biomol/seq-id-map``), fetched per the README.

Destinations follow the same ``DATA_ROOT`` / ``OUTPUT_ROOT`` layout the db configs read.
"""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

CCD_URL = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
SABDAB_URL = "https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/summary/all/"
# wwPDB mmCIF rsync mirror (divided into 2-char subdirs). Full PDB mmCIF is ~90 GB+.
MMCIF_RSYNC = "rsync.rcsb.org::ftp_data/structures/divided/mmCIF/"
OPENFOLD_PORTAL = "https://portal.openfold.omsf.io/"


def _curl(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest`` with curl (resumable, follows redirects)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(  # noqa: S603
        ["curl", "-fL", "--retry", "3", "-C", "-", "-o", str(dest), url],  # noqa: S607
        check=True,
    )


def download_ccd(output_root: Path) -> Path:
    """Fetch the wwPDB CCD and split it into one ``<COMP_ID>.cif`` per component.

    The ccd build is file-per-component (``OUTPUT_ROOT/materials/raw/ccd/components``),
    so the single concatenated ``components.cif`` is split on its ``data_`` blocks.
    Returns the components directory.
    """
    raw_dir = output_root / "materials" / "raw" / "ccd"
    gz_path = raw_dir / "components.cif.gz"
    components = raw_dir / "components"
    _curl(CCD_URL, gz_path)
    components.mkdir(parents=True, exist_ok=True)

    n = 0
    out = None
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("data_"):
                if out is not None:
                    out.close()
                comp_id = line.strip()[len("data_"):]
                out = (components / f"{comp_id}.cif").open("w", encoding="utf-8")
                n += 1
            if out is not None:
                out.write(line)
    if out is not None:
        out.close()
    return components


def download_sabdab(data_root: Path) -> Path:
    """Fetch the SabDab antibody summary TSV (used by seq_cluster). Returns its path."""
    dest = data_root / "external" / "SabDab" / "sabdab_summary_all.tsv"
    _curl(SABDAB_URL, dest)
    return dest


def download_mmcif(data_root: Path, *, confirmed: bool) -> Path:
    """Rsync the full wwPDB mmCIF mirror (~90 GB+). Requires ``confirmed=True``.

    Mirrors the divided (2-char subdir) layout; the pdb/cif build reads
    ``DATA_ROOT/mmcif_files_latest/mmcif_files``.
    """
    dest = data_root / "mmcif_files_latest" / "mmcif_files"
    if not confirmed:
        msg = (f"mmCIF is a ~90 GB+ rsync from {MMCIF_RSYNC}. "
               f"Re-run with --yes to start it (target: {dest}).")
        raise RuntimeError(msg)
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(  # noqa: S603
        ["rsync", "-rlptz", "--delete", MMCIF_RSYNC, str(dest)],  # noqa: S607
        check=True,
    )
    return dest


def download_openfold(data_root: Path, *, confirmed: bool) -> Path:  # noqa: ARG001
    """Document the OpenFold3 distillation download -- it is TB-scale and portal-hosted.

    The distillation sets are served from an interactive portal, not a single bulk URL,
    and are enormous, so this deliberately does not auto-fetch: it prints where to get
    them and the destination layout the db/distillation configs expect. Even with
    ``confirmed=True`` it only prints instructions.
    """
    dest = data_root / "openfold_distillation"
    msg = (
        f"OpenFold3 distillation sets are TB-scale and served from an interactive "
        f"portal, so they are not auto-downloaded.\n"
        f"  1. Get the monomer (long/short), RNA, and disordered sets from:\n"
        f"       {OPENFOLD_PORTAL}\n"
        f"  2. Place them under: {dest}\n"
        f"     (monomer_distillation_sets_v2/, rna_distillation_set/, disordered_set/, "
        f"fasta/ -- see db/distillation/*.yaml for exact subpaths).\n"
        f"  Only build these if you truly need to rebuild them (the project rule is to "
        f"validate small vs the existing DBs, not full-build)."
    )
    raise RuntimeError(msg)
