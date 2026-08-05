# Manual mmCIF fixes (pdb/cif)

A set of **53 PDB entries** cannot be ingested from the *current* wwPDB mmCIF snapshot:
they either error out during the `pdb/cif` build or would build wrongly. The original
BioMol build handled them by **substituting an older, known-good mmCIF** for each entry
before ingest. Reproducing production decode-identically requires the same substitution,
so it is a first-class, documented step here — not a silent workaround.

- **The id list** lives in [`db/pdb/manual_cif_fixes.txt`](../db/pdb/manual_cif_fixes.txt) (53 entries).
- **Apply it** with `structcooker fix-cif --source <corrected-cif-dir>` before
  `structcooker build pdb/cif`.
- **`structcooker inspect`** reports how many substitutions are in place
  (`[manual cif fixes] PENDING n/53 …`).

## Why these entries need it

Verified end-to-end on `1ai0` (a 10-model NMR ensemble):

- The same non-polymer ligand (`IPH`, phenol) is **re-numbered per model** in
  `_atom_site`: `auth_seq_id = 22` in models 1,2,3,6,8 and `31` in models 4,5,7,9,10.
- `_pdbx_nonpoly_scheme` records the ligand **once** (`pdb_seq_num = auth_seq_num = 22`).
- The build matches each atom to its scheme residue by `auth_seq_id`. For the models
  numbered 31, the atom's `auth_seq_id` is absent from the scheme (which only has 22), so
  `_scatter_atom_site_coords` raises `Auth_idx 31 not found in scheme … residue IPH`.

The failure is intrinsic to the mmCIF file: the per-model renumbering lives in
`_atom_site`, the single number in `_pdbx_nonpoly_scheme` — both from the file, neither
from the CCD. Most of the 53 share this NMR-ligand pattern.

### What it is *not*

- **Not a rewrite/port regression.** The scheme parser, atom `auth_idx` construction, and
  scatter/merge logic are byte-identical from the repo's first cif commit (2026-03) through
  the pre-refactor legacy (2026-04) to today. Every version fails on these files the same way.
- **Not a CCD difference.** In the failing check the CCD is consulted only to test whether
  the residue is named `WATER` (the sole tolerated skip). `IPH` is never `WATER`, so no CCD
  build changes the outcome; the `22`-vs-`31` mismatch is pure mmCIF.

The original build simply never fed these files to the parser: it replaced them first.

## How production did it

The legacy pipeline's `scripts/manually_fix_cif.py` copied each error entry's cif from an
older snapshot (`BioMolDB_2024Oct21/cif/cif_raw/`, divided `<id[1:3]>/` layout) over the
build's raw cif. The error-item list was `…/cif/error_items/manually_fixed/` — the 53 ids
captured here. All entries that fail the current clean build are a subset of this list.

## Reproducing it

1. **Provide the corrected cifs.** The older snapshot is a provided external input (not
   redistributed here). Point `--source` at a directory holding a good cif per id, either
   divided (`<id[1:3]>/<id>.cif[.gz]`) or flat (`<id>.cif[.gz]`).
2. **Apply before building:**
   ```bash
   structcooker fix-cif --source /path/to/corrected-cif-snapshot   # overlays into the mmCIF dir
   structcooker inspect                                            # [manual cif fixes] OK
   structcooker build pdb/cif
   ```
   `fix-cif` decompresses `.gz` sources to `<id>.cif`, writes an `.manual_cif_fixes_applied`
   marker so `inspect` can see the state, and `--dry-run` reports without writing.

> **Note on `DATA_ROOT`.** `fix-cif` writes into `DATA_ROOT/mmcif_files_latest/mmcif_files`.
> Run reproduction against a `DATA_ROOT` you own; do not overlay a shared/read-only mmCIF
> mirror in place.

## Current build state (decision: build without the substitution for now)

The corrected snapshot (`BioMolDB_2024Oct21`) is **not present on this cluster**, and the
error-marker copies that are (`…/error_items/manually_fixed/`) are the *broken* files, not
fixes — a 25/25 isolated rebuild of them fails identically. So the accepted clean build
(`cif_pdb.lmdb`, **233,579** entries) was produced from the current wwPDB mmCIF **without**
the substitution:

- **25** of the 53 error out and are absent (of these, 9 would reach the training views).
- **28** of the 53 build from the *current* mmCIF and are present, but they were **not** run
  through production's substitution, so those 28 are not guaranteed byte-identical to
  production. (Everything outside the 53 is decode-identical — verified on shared keys.)

`structcooker fix-cif` is wired and ready: once the corrected snapshot is obtained, apply it
and rebuild to fold all 53 in production-faithfully.

## Alternative: fix the parser instead of substituting

Rather than substitute files, the entries could be built from the *current* mmCIF by making
the atom→scheme match **model-invariant** — union the per-model `auth_seq` numbers into the
scheme, or match by residue identity rather than raw `auth_seq_id`. This is self-contained
(no external snapshot) but would **not** reproduce production's exact bytes for these entries
(different mmCIF revision). Tracked as future work.
