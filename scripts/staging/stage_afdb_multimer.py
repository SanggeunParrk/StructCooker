"""Stage the AFDB multimer structures into materials/raw, the tree the AFM builds read.

Two sources, handled differently because of who owns them and what they contain:

* Homodimers -- the high-confidence set from the AFDB complex paper (Han et al.), shared by
  the first author as ``high_conf_homodimers_manuscript_v1.tar.gz`` and extracted under a
  personal scratch path. 1,750,755 PDBs. These are MOVED (same filesystem, so each move is
  a rename); the tarball stays where it is as the source of truth.
* Heterodimers -- the 80,248 ``passes_quality_threshold`` entries fetched one by one by
  hatapakacha5 into ``AFDB_heterodimer/cif``. That directory is someone else's work and not
  writable to us, so these are COPIED. Their ``.cif`` files are gzip streams (the server sent
  gzip and the fetch saved it verbatim), and the CIF reader picks its decoder by extension, so
  the copies are named ``.cif.gz``. The bytes are unchanged.

Names are normalised: 35,254 homodimer files are ``AF_<n>`` rather than ``AF-<n>`` -- a quirk
of the release (no id occurs in both forms). The entry key comes from the filename and the
metadata uses ``AF-`` throughout, so a mixed tree would silently fail to join. The original
name is kept in SOURCE.tsv.

Both trees shard on the last three digits of the entity number, like teddymer's DimerIndex
and the MSA tree's seq_id shards, so a path is computable from an id.

Resumable: a move that already happened is simply absent from the source directory, and a
copy whose destination already has the right size is skipped.

Usage:
    python scripts/staging/stage_afdb_multimer.py --out /data/shared/cssb_data/BioMol/materials/raw/afdb_multimer
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
from pathlib import Path

HOMO_SRC = Path(
    "/data/psk6950/external_source/AFDB/multimer/test/high_conf_homodimers_manuscript_v1",
)
HET_SRC = Path("/data/shared/cssb_data/AFDB_heterodimer/cif")
_ID = re.compile(r"^AF[-_](\d+)-model_v1\.(pdb|cif)$")


def _dest(root: Path, number: str, name: str) -> Path:
    return root / number[-3:] / name


def stage_homodimers(out: Path) -> int:
    """Move every homodimer PDB into ``out/homodimer/pdb``; return how many are staged."""
    root = out / "homodimer" / "pdb"
    for shard in range(1000):
        (root / f"{shard:03d}").mkdir(parents=True, exist_ok=True)
    log = (out / "homodimer" / "SOURCE.tsv.partial").open("a")
    with os.scandir(HOMO_SRC) as entries:
        for entry in entries:
            m = _ID.match(entry.name)
            if m is None or m.group(2) != "pdb":
                continue
            number = m.group(1)
            name = f"AF-{number}-model_v1.pdb"
            dest = _dest(root, number, name)
            size = entry.stat().st_size
            # Record before renaming: after the rename the original name is gone.
            log.write(f"AF-{number}\t{entry.name}\t{entry.path}\t{size}\n")
            log.flush()
            Path(entry.path).rename(dest)
    log.close()
    return _finish_source(out / "homodimer", "entity_id\toriginal_name\toriginal_path\tbytes")


def stage_heterodimers(out: Path) -> int:
    """Copy every heterodimer mmCIF into ``out/heterodimer/cif`` as .cif.gz."""
    root = out / "heterodimer" / "cif"
    for shard in range(1000):
        (root / f"{shard:03d}").mkdir(parents=True, exist_ok=True)
    rows = []
    with os.scandir(HET_SRC) as entries:
        for entry in entries:
            m = _ID.match(entry.name)
            if m is None or m.group(2) != "cif":
                continue
            number = m.group(1)
            dest = _dest(root, number, f"AF-{number}-model_v1.cif.gz")
            size = entry.stat().st_size
            if not dest.exists() or dest.stat().st_size != size:
                shutil.copy2(entry.path, dest)
            rows.append(f"AF-{number}\t{entry.name}\t{entry.path}\t{size}\n")
    with (out / "heterodimer" / "SOURCE.tsv").open("w") as fh:
        fh.write("entity_id\toriginal_name\toriginal_path\tbytes\n")
        fh.writelines(sorted(rows))
    return len(rows)


def _finish_source(part_dir: Path, header: str) -> int:
    """De-duplicate the append-only move log (a resumed run re-logs nothing already moved)."""
    partial = part_dir / "SOURCE.tsv.partial"
    rows = sorted(set(partial.read_text().splitlines())) if partial.exists() else []
    (part_dir / "SOURCE.tsv").write_text(header + "\n" + "".join(r + "\n" for r in rows))
    partial.unlink(missing_ok=True)
    return len(rows)


def write_readme(out: Path, n_homo: int, n_het: int) -> None:
    (out / "README.md").write_text(
        f"""# afdb_multimer — raw input for the AFM (AFDB multimer) DB

High-confidence predicted dimers from the AFDB complex expansion (Han et al., "AlphaFold
Database expands to proteome-scale quaternary structures", 2026), staged by
`scripts/staging/stage_afdb_multimer.py`.

| part | entries | format | origin |
|---|---:|---|---|
| `homodimer/pdb/` | {n_homo:,} | PDB | first author's `high_conf_homodimers_manuscript_v1.tar.gz` |
| `heterodimer/cif/` | {n_het:,} | mmCIF, gzip | EBI per-entry fetch of `passes_quality_threshold=true` |

## What "high confidence" means here

Homodimers: the manuscript set. The released `homodimer_metadata.csv` carries no quality
flag and no pLDDT column, so the paper's selection (ipSAE_min >= 0.6 plus a pLDDT/clash
filter) cannot be re-applied from metadata -- ipSAE alone gives 2,593,499. The manuscript
bundle (1,750,755; paper reports 1,754,242) is the set, and it sits inside ipSAE_min >= 0.6
with only 7 exceptions.

Heterodimers: the recalibrated AFDB release (80,248), not the paper's 56,959 "tentatively
high-confidence" -- the paper deferred heterodimer calibration to a later release.

## Layout

    homodimer/pdb/<last 3 digits>/AF-<n>-model_v1.pdb
    heterodimer/cif/<last 3 digits>/AF-<n>-model_v1.cif.gz
    */SOURCE.tsv   entity_id, original name, original path, bytes

Names are normalised to `AF-<n>`: the release spells 35,254 homodimers `AF_<n>`. Heterodimer
files were `.cif` holding gzip bytes; they are `.cif.gz` here so readers that dispatch on
extension decode them. Contents are byte-identical to the sources.

B-factor holds pLDDT. Homodimer chains are identical sequences; heterodimer chains differ.
""",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    n_het = stage_heterodimers(args.out)
    print(f"heterodimer copied {n_het}", flush=True)
    n_homo = stage_homodimers(args.out)
    print(f"homodimer   moved  {n_homo}", flush=True)
    write_readme(args.out, n_homo, n_het)

    left = sum(1 for e in os.scandir(HOMO_SRC) if e.name.endswith(".pdb"))
    if left:
        print(f"WARNING: {left} homodimer PDBs still in the source directory")


if __name__ == "__main__":
    main()
