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

Units (`unit/` stages of the DAG, [docs/unit-db.md](../docs/unit-db.md)): `lmdb/unit/pdb`, built 2026-11-10 from `lmdb/pdb/cif/cif_pdb.lmdb`:
944,570 units from 247,987 of the 249,676 PDB entries (the other 1,689 give no unit): 850,925 contact units (the same as the
earlier contact-only `lmdb/interface/pdb`), 89,904 monomer units (85,898 protein, 4,006 NA) and 3,741 base-paired NA duplexes
without a contact edge. Split by release date: train 547,053 (167,321 entries) · valid 92,091 (16,429) · test 262,839 (47,011) ·
none 42,587 (resolution rule: no resolution 16,392, > 4.5 A 18,408, > 9 A 7,787; no unit lacks a date). 160.6 GB. Built stage by stage
(`structcooker build unit/<stage> --depends-on`), not by `build-all`, which refuses the existing roots. `lmdb/interface/pdb` is the superseded
contact-only build.
`scripts/unit/validate_db.py` on 3,000 random units: every check passes except `residue_beyond_sequence` on `9jjb_1_1_0` (two
protein nodes keyed with the 94-residue sequence of chains B / D, while the node's residues are the 552 resolved residues of the
605-residue chains, cif_idx up to 597); the same unit with the same numbers is in the earlier build, so monomers did not cause it. Cause
not found yet; in that entry `auth_asym_id` (A, C, B, D) is permuted against `chain_id` (A_1, B_1, C_1, D_1), a lead for node-to-residue pairing.

Derived sets: train_20210930 167,912 · train_20260301 233,579 · mpnn 484,383 · valid1 21,527 ·
valid2 617 · cif_chain 4,162,929 chains · CCD 50,782 components (+ inference CCD 50,781).

## Known gaps

- **Units**: the AF3 crystallisation-aid list (16 CCD ids) was transcribed from memory of the AF3
  SI and is not checked against the PDF; detergents / lipids are excluded by our choice; no buried SASA, no
  experimental method column.
- **PDB templates** are the 2026-07-30 build: the 293 entries (51,895 chains) `cif_chain`
  gained on 2026-09-30 are not template candidates. Templates of every set are restricted
  to hits released by 2021-09-30.
- **Failed inputs** (`db_configs/release-exceptions.json`): 5 PDB entries the parser cannot map, and
  840 OpenFold `alignment.npz` that are empty or corrupt in the release (long 667, short 173).
- **OFD long**: 80 entries have no model in the OpenFold release.
- **pLDDT**: none published for the OFD RNA and disordered sets.

## Checks (2026-10-13)

- `structcooker release-check`: every stage present, self-contained and consistent.
- `structcooker e2e-sample`: a sample of every source rebuilds, through all 95 stages,
  record-for-record identical to this release for every structure, MSA and CCD record.

## Recent changes

- 2026-11-10: the interface DB is renamed the unit DB (`unit/` stages, `lmdb/unit/pdb`, `instructions/unit`) and gains
  monomer units (polymer chains in no contact unit, `monomer` in `params/unit.yaml`); `lmdb/unit/pdb` built the same day (944,570 units).
- 2026-10-21: interface-unit DB of the PDB (`lmdb/interface/pdb`): units of 2..N touching molecules (ions and ligands
  as nodes, antibody / TCR pairs and base-paired duplexes always bundled, waters stored), date split
  train / valid / test; since 2026-11 built by the `unit/` stages of the DAG (formerly `structcooker interface`).
- 2026-10-13: one root (`BIOMOL_ROOT`) holding raw, reference and intermediate inputs; teddymer
  MSAs and the seq_id seed became supplied inputs (no search in the release build); OFD long
  shards gathered inside their DBs; PDB templates and FoldBench assets moved in; cluster
  settings from the environment (`scripts/cluster/cssb.env`); one build surface.
- 2026-10-12: OFD long MSAs gained the 14 entries recovered 2026-09-14; CCD rebuilt with the
  current code (only 5j2p changed, patched downstream).
