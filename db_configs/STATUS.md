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

Units (`unit/` stages of the DAG, [docs/unit-db.md](../docs/unit-db.md)): `lmdb/unit/pdb`, built 2026-11-10 from `lmdb/pdb/cif/cif_pdb.lmdb` and
patched 2026-11-11 (below): 944,828 units from 247,987 of the 249,676 PDB entries (the other 1,689 give no unit): 851,207 contact units,
89,871 monomer units (85,869 protein, 4,002 NA) and 3,750 base-paired NA duplexes without a contact edge. Split by release date: train
547,194 (167,321 entries) · valid 92,105 (16,429) · test 262,935 (47,011) · none 42,594 (resolution rule: no resolution 16,399, > 4.5 A
18,408, > 9 A 7,787; no unit lacks a date). Built stage by stage (`structcooker build unit/<stage> --depends-on`), not by `build-all`, which
refuses the existing roots. `lmdb/interface/pdb` is the superseded contact-only build.
Patch of 2026-11-11 (a node-merge bug): a ligand covalently bonded to two polymer chains used to fuse both chains into one node, keyed by the
sequence of one of them (`9jjb_1_1_0`: a 605- and a 94-residue chain in one node keyed with the 94-residue sequence, found by
`validate_db.py` as `residue_beyond_sequence`). 228 of the 249,676 entries had such a ligand in an assembly the DB reads
(`scripts/unit/bridge_entries_scan.py`); their 510 stored units became 768 (`structcooker patch unit/unit_pdb --keys`, then the 44 old ids that
no rebuilt unit overwrote were deleted with `scripts/unit/drop_unit_ids.py`; manifest, manifest_full and split rebuilt; the pre-patch tables are in
`lmdb/unit/pdb/.before_bridge_fix/`). A scan of all 944,828 units (`scripts/unit/merged_node_scan.py`) finds no node with more than one chain, and
`scripts/unit/validate_db.py` on 3,000 random units passes every check. `release-check` (2026-11-12) passes `unit/*` (944,828 units, 200 sampled); `e2e-sample` does not sample the unit stages yet.

Re-Glyco (`glycan/` stages, [docs/reglyco-loader.md](../docs/reglyco-loader.md), [docs/reglyco-compact.md](../docs/reglyco-compact.md)): `lmdb/glycan/reglyco/reglyco.lmdb`,
built 2026-11-11/12 from the RFD4-Proteina Re-Glyco component (Hub commit `aa86df1e`, sha256-verified): 2,992,448 records of 11,160 PDB entries as
BioMol records stored as deltas (codec `biomol_zstd_delta`, schema K; read with `DeltaReader`), 14.95 GB (the shipped component is 1,574.6 GB). The Hub's per-record
columns are in `metadata/reglyco_examples.parquet`. Built stage by stage like the unit DB (`structcooker build glycan/<stage> --depends-on`, with
`DATACOOKER_ADMISSION_START=32`). Checked: every record decodes (2,992,448 / 2,992,448); 5,000 random records equal a fresh conversion of their CIF (5,000 / 5,000);
a 3,000-group sweep converted without failure. `release-check` (2026-11-12) passes `glycan/*` (2,992,448 records, 200 sampled); `e2e-sample` does not sample the Re-Glyco stages yet.

RFD4-Proteina sets (`rfd4/` stages, [docs/rfd4-sets.md](../docs/rfd4-sets.md)): `lmdb/rfd4/<set>/<set>.lmdb`, built 2026-11-12 from Hub commit `aa86df1e` as ordinary BioMol records (schema M, the PDB's own
recipe): `omol` 53,656,064 (of 53,664,160: 8,096 name a component no CCD or record defines, or are empty; `metadata/rfd4_omol_unconvertible.tsv`), `rna_monomer` 371,017,
`synthetic_ncaa` 1,463,222, `ncaa_scaffolds` 288,390, `antibody_distill` 75,576, `na_complex` 20,261, `tcr_distill` 3,941 (the other six lose nothing). The Hub's per-record columns are in
`metadata/rfd4_<set>_examples.parquet`. Checked: `release-check` passes every `rfd4/*` stage; a sample of each decodes with residue counts equal to the Hub's; the sampled end-to-end run
(`e2e-sample ... --manifest db_configs/MANIFEST_rfd4.yaml`, 40 records per set) rebuilds every sampled record identically to the release.
PTM / glycan / cyclic tags of PDB chains: `metadata/pdb_chain_tags.csv` (171,952 chain rows from the scan of `cif_pdb`; 32 entries the scan could not read are missing), a separate table, not in any DB.
`release-check` (2026-11-12, 127 stages) reports 10 stages with problems, none of them above: `design/monomer_{long,short}` and their `_clusters/_items/_stats` (the design set is not built), `metadata/ab_seq`, `metadata/cdr_region`.

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

- 2026-11-12: RFD4-Proteina structure sets (`rfd4/`: omol, rna_monomer, na_complex, antibody_distill, tcr_distill, synthetic_ncaa, ncaa_scaffolds) built as BioMol records; the CCD step accepts a component the record itself defines; `metadata/pdb_chain_tags.csv`.

- 2026-11-10: the interface DB is renamed the unit DB (`unit/` stages, `lmdb/unit/pdb`, `instructions/unit`) and gains
  monomer units (polymer chains in no contact unit, `monomer` in `params/unit.yaml`); `lmdb/unit/pdb` built the same day (944,570 units; 944,828 after the patch of 2026-11-11).
- 2026-10-21: interface-unit DB of the PDB (`lmdb/interface/pdb`): units of 2..N touching molecules (ions and ligands
  as nodes, antibody / TCR pairs and base-paired duplexes always bundled, waters stored), date split
  train / valid / test; since 2026-11 built by the `unit/` stages of the DAG (formerly `structcooker interface`).
- 2026-10-13: one root (`BIOMOL_ROOT`) holding raw, reference and intermediate inputs; teddymer
  MSAs and the seq_id seed became supplied inputs (no search in the release build); OFD long
  shards gathered inside their DBs; PDB templates and FoldBench assets moved in; cluster
  settings from the environment (`scripts/cluster/cssb.env`); one build surface.
- 2026-10-12: OFD long MSAs gained the 14 entries recovered 2026-09-14; CCD rebuilt with the
  current code (only 5j2p changed, patched downstream).
