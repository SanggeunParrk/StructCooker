# Release status

The current BioMol release (`BIOMOL_ROOT`, default `/data/shared/cssb_data/BioMol_clean`).
How it is built: [docs/e2e-build.md](../docs/e2e-build.md). The tree's own catalogue, sizes and
check results: `$BIOMOL_ROOT/docs/` (`structcooker release-docs`). Earlier notes:
[docs/history/status-ledger-2026-10.md](../docs/history/status-ledger-2026-10.md).

## Contents

MSA and template DBs are keyed by seq_id, so they count distinct sequences, not structures.

| set | structures (cif / cif_attached) | MSAs (+ depth caps) | templates | b_factor |
|---|---|---|---|---|
| PDB | 249,676 / 249,644 | 178,249 (+d16k/d2k/d512); RNA 6,572 (+d2k/d512) | 1,007,928 chains; 176,579 seq_ids | experimental |
| OFD short | 430,418 / 430,418 | 430,245 (+d2k/d512/d5k) | 430,370 | pLDDT |
| OFD long | 16,099,404 / 16,099,404 | 16,098,808 (+d2k/d512) | 16,080,773 | pLDDT |
| OFD rna | 126,778 / 126,778 | 126,751 (+d2k/d512/d5k) | — | 0 (not published) |
| OFD disordered | 28,567 / 28,567 (+contacts) | 19,649 (+d2k/d512) | 78,176 | 0 (not published) |
| Teddymer | 510,454 / 510,454 | 999,853 (+d2k) | 998,878 | pLDDT |
| AFDB homodimer | 1,750,755 / 1,750,755 | 1,721,635 (+d2k), shared with heterodimers | 1,720,477, shared | pLDDT |
| AFDB heterodimer | 80,248 / 80,248 | (above) | (above) | pLDDT |

Derived sets: train_20210930 167,912 · train_20260301 233,579 · mpnn 484,383 · valid1 21,527 ·
valid2 617 · cif_chain 4,162,929 chains · CCD 50,782 components (+ inference CCD 50,781).

## Known gaps

- **PDB templates** are the 2026-07-30 build: the 293 entries (51,895 chains) `cif_chain`
  gained on 2026-09-30 are not template candidates. Templates of every set are restricted
  to hits released by 2021-09-30.
- **Failed inputs** (`db/release-exceptions.json`): 5 PDB entries the parser cannot map, and
  840 OpenFold `alignment.npz` that are empty or corrupt in the release (long 667, short 173).
- **OFD long**: 80 entries have no model in the OpenFold release.
- **pLDDT**: none published for the OFD RNA and disordered sets.

## Checks (2026-10-13)

- `structcooker release-check`: every stage present, self-contained and consistent.
- `structcooker e2e-sample`: a sample of every source rebuilds, through all 94 stages,
  record-for-record identical to this release for every structure, MSA and CCD record.

## Recent changes

- 2026-10-13: one root (`BIOMOL_ROOT`) holding raw, reference and intermediate inputs; teddymer
  MSAs and the seq_id seed became supplied inputs (no search in the release build); OFD long
  shards gathered inside their DBs; PDB templates and FoldBench assets moved in; cluster
  settings from the environment (`scripts/cluster/cssb.env`); one build surface.
- 2026-10-12: OFD long MSAs gained the 14 entries recovered 2026-09-14; CCD rebuilt with the
  current code (only 5j2p changed, patched downstream).
