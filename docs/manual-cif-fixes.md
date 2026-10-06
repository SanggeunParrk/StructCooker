# Historical mmCIF substitutions

An earlier BioMol snapshot used 53 substitutions listed in
[`manual_cif_fixes.txt`](../db/pdb/manual_cif_fixes.txt). This is historical provenance,
not a required preprocessing step for every PDB release. The September 10 recovery
used the provided snapshot directly; its five rejected records and causes are in the
[recovery report](recovery-2026-09-10.md).

When reproducing the older snapshot, obtain its corrected files and apply them to a
copy of the raw input that you own:

```bash
export MMCIF_ROOT=/path/to/your/copied/mmcif
structcooker fix-cif --source /path/to/corrected-snapshot --dry-run
structcooker fix-cif --source /path/to/corrected-snapshot
```

The command reads flat or divided sources, compressed or uncompressed. It replaces
the existing flat/divided `<id>.cif.gz` consumed by the CIF recipe, staging each file
before an atomic rename. An invalid compressed source leaves the original intact.
Duplicate compressed destinations are rejected. The marker
`.manual_cif_fixes_applied` records the union of applied IDs; it is bookkeeping, not
a content checksum or a scientific validation. General preflight no longer treats
these historical substitutions as pending work for a different snapshot.

Never apply replacements to the shared production snapshot. Without `MMCIF_ROOT`,
the destination is `DATA_ROOT/BioMol/materials/raw/cif`.

## Not needed with the current parser (2026-10-01)

All 53 listed records were compared with production `BioMol/lmdb/pdb/cif/cif_pdb.lmdb`. The 48
that build from the current snapshot are identical to production in every array; the other 5
(2g10, 2icy, 2q44, 4xq2, 9gdy) are absent from production as well. The substitutions are
therefore not part of the BioMol_clean build; `fix-cif` stays for reproducing the older snapshot.
