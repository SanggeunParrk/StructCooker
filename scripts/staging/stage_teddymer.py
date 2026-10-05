"""Copy the Teddymer selection into materials/raw, the tree db/teddymer/cif builds from.

Teddymer ships 10,089,503 dimer PDBs in two flat directories under a personal scratch
path. We build 510,454 of them (see db/MANIFEST_teddymer.yaml for why), so this copies
exactly that set into the shared materials tree alongside the other raw inputs -- real
files, so the build no longer depends on the scratch copy surviving.

Selection = the source paper's filter over the non-singleton cluster representatives,
which are the only entries with published interface metrics:
    InterfaceLength > 10,  AvgIntPAE < 10,  AvgIntPlddt > 70

Usage:
    python scripts/staging/stage_teddymer.py \
        --teddymer-root /data/psk6950/external_source/teddymer \
        --out /data/shared/cssb_data/BioMol/materials/raw/teddymer
"""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

# Interface-quality thresholds. Reproducing these over nonsingletonrep_metadata.tsv gives
# 510,454 representatives, matching the source paper exactly [verified 2026-09-22].
MIN_INTERFACE_LENGTH = 10.0
MAX_INTERFACE_PAE = 10.0
MIN_INTERFACE_PLDDT = 70.0

# Shard on the key, the way materials/raw/cif shards on the PDB id's middle two characters:
# the DimerIndex's last three digits give 1,000 directories of ~510 files, so a path is
# computable from an entry id instead of needing a scan. scan_paths recurses into them.
SHARD_WIDTH = 3

_DIMER_SUBDIRS = ("teddymer_dimer_pdbs_1", "teddymer_dimer_pdbs_2")


def passing_dimer_indices(metadata_path: Path) -> set[str]:
    """DimerIndex of every representative clearing the interface filter."""
    passing: set[str] = set()
    with metadata_path.open() as handle:
        next(handle)  # header
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            if (
                float(parts[4]) > MIN_INTERFACE_LENGTH
                and float(parts[5]) < MAX_INTERFACE_PAE
                and float(parts[6]) > MIN_INTERFACE_PLDDT
            ):
                passing.add(parts[0])
    return passing


def stage(teddymer_root: Path, out_dir: Path) -> tuple[int, int]:
    """Copy every selected PDB into ``out_dir``; return (copied, selected)."""
    selected = passing_dimer_indices(
        teddymer_root / "teddymer_data" / "teddymer" / "nonsingletonrep_metadata.tsv",
    )
    pdb_dir = out_dir / "pdb"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    for shard in range(10**SHARD_WIDTH):
        (pdb_dir / f"{shard:0{SHARD_WIDTH}d}").mkdir(exist_ok=True)
    _write_readme(out_dir, teddymer_root, len(selected))

    record = (out_dir / "SOURCE.tsv").open("w")
    record.write("dimer_index\tsource_path\tbytes\n")
    copied = 0
    for sub in _DIMER_SUBDIRS:
        source_dir = teddymer_root / "dimer_pdbs" / sub
        # scandir streams: a 5M-entry directory must not be materialised as a list.
        with os.scandir(source_dir) as entries:
            for entry in entries:
                if not entry.name.endswith(".pdb"):
                    continue
                # The filename carries the TED pair suffix the cluster ids do not;
                # the DimerIndex prefix is the key all three files share.
                index = entry.name.split("DI_", 1)[0]
                if index not in selected:
                    continue
                dest = pdb_dir / index.zfill(SHARD_WIDTH)[-SHARD_WIDTH:] / entry.name
                size = entry.stat().st_size
                # Resumable: an already-complete copy is left alone, so a killed run
                # just continues. A short file is re-copied rather than trusted.
                if not dest.exists() or dest.stat().st_size != size:
                    shutil.copy2(entry.path, dest)
                record.write(f"{index}\t{entry.path}\t{size}\n")
                copied += 1
    record.close()
    return copied, len(selected)


def _write_readme(out_dir: Path, teddymer_root: Path, selected: int) -> None:
    """Leave the folder able to explain itself without this script."""
    (out_dir / "README.md").write_text(
        f"""# teddymer — raw input for db/teddymer/cif

{selected:,} TED domain-pair structures, copied by
`scripts/staging/stage_teddymer.py` from `{teddymer_root}`.

## What these are

Each file is TWO TED domains of ONE AFDB v4 model, split into chains A/B -- a
domain-domain interface, not a protein-protein complex. A domain can be
sequence-discontinuous, so a chain's residue numbers may jump; the numbering is the
original model's and a jump marks a domain boundary, not a gap in the structure.
ATOM records only; the B-factor column holds pLDDT.

## Why only {selected:,} of 10,089,503

Teddymer clusters into 3,556,223 clusters (587,687 non-singleton + 2,968,536
singleton). Interface quality is published only for non-singleton cluster
representatives, so only those can be filtered. The source paper's filter --
InterfaceLength > {MIN_INTERFACE_LENGTH:.0f}, AvgIntPAE < {MAX_INTERFACE_PAE:.0f},
AvgIntPlddt > {MIN_INTERFACE_PLDDT:.0f} -- leaves {selected:,} of them.
See db/MANIFEST_teddymer.yaml for the sets that were measured and rejected.

## Layout

    pdb/<last 3 digits of DimerIndex>/<DimerIndex>DI_<accession>_<TED pair>.pdb
    SOURCE.tsv   dimer_index, source path, bytes -- provenance and integrity

The DimerIndex prefix is the join key against teddymer's `cluster.tsv` and
`nonsingletonrep_metadata.tsv`; neither of those carries the TED suffix.
""",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--teddymer-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    copied, selected = stage(args.teddymer_root, args.out)
    print(f"selected {selected}  copied {copied}")
    if copied != selected:
        # A representative with no PDB on disk is a staging bug, not a filter outcome.
        print(f"WARNING: {selected - copied} selected entries had no source file")


if __name__ == "__main__":
    main()
