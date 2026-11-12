# Re-Glyco, stored compactly

The RFD4-Proteina Re-Glyco set (`train_datasets/streaming-latest/synthetic/glycans/2026_09_27_reglyco`, see
[raw materials](raw-materials.md)) is GlycoShape glycans rebuilt on PDB glycosylation sites: **2,992,448 records** of
11,160 PDB entries, **1,573.6 GB** as shipped (10.85 TB decompressed). This page records why that is mostly repeated
protein, how small a random-access copy can be, and the layout to build. Everything is measured on one shipped shard
(`rs_00005.bin`, 700 records of 5 entries) and the Hub's two tables, `blob/index.parquet` and `examples.parquet`.

**The release has one value format, BioMol** (`pdb/cif` schema A: `{assembly_dict, metadata_dict}`), so the compact copy
has to hold BioMol records, not CIF text. Measured on BioMol records:

- Our `pdb/cif` recipe converts all 700 records (0 failures) once a small adapter derives the tables the annotated CIF
  lacks ([below](#as-biomol-records)); the conversion keeps every atom with coordinates (checked on 1,059 assemblies).
- One BioMol value per record (`to_bytes`, zstd 6) is **0.75 of the shipped CIF: about 1,180 GB for the set (900 GB for the
  training view)**. Converting alone does not make this small; the records repeat the protein.
- A **delta against a base record of the same entry** stores the set in **about 7.5 GB (5.8 GB for the training view)** at
  zstd level 19, **12 GB at level 9**, with every record decoding to a BioMol dict equal to the converted one. Level 9 is
  not enough here: unlike CIF text, BioMol arrays shift when the glycan changes, and only level 19 finds the match.
- The earlier measurement on raw CIF text ([below](#layouts-on-raw-cif-text-the-earlier-measurement)) gave 8.4 GB at level 9
  and 7.0 GB at level 19. The BioMol form costs nothing extra at level 19, but the build is about 100 times more expensive
  (the CIF has to be parsed and run through the recipe: about 1,100 core-hours, not 10).

Code: `scripts/analysis/reglyco_biomol.py` (+ `.sbatch`), `reglyco_biomol_extrapolate.py`, and for the CIF text variant
`reglyco_compact.py`, `reglyco_extrapolate.py`. Outputs: `$BIOMOL_ROOT/materials/intermediate/reglyco_biomol/rs_00005/`
and `.../reglyco_compact/`.

## What differs between records

A record is one annotated mmCIF named `<xx>/<pdb>/<pdb>_<site>_<glycan>_conf<NNN>.cif`: one PDB entry, one glycosylation
site, one GlycoShape glycan, one conformer. The index gives, for every record, `shard`, `offset`, `length`,
`uncompressed_length` and the sha256 of the decompressed CIF.

- The whole protein context is rewritten into every record: 14,376 lines for `2pmv_A293_GS00178_conf000`.
- Two conformers of one site differ in **83 of those lines**, all glycan `HETATM` rows (coordinates only).
- Another site or glycan of the same entry replaces the glycan rows.
- The population: 11,160 entries, 499,626 (entry, site, glycan) groups, 6.0 conformers per group, 268 records per entry
  (up to 3,261). 11,064 entries lie entirely in one shard. Records are 3.8 MB at the median, 9.4 MB at the 99th percentile
  and 13.7 MB at most.

## As BioMol records

`reglyco_biomol.py convert` runs the `pdb/cif` recipe (`workflows/ingest/cif.py`, targets `assembly_dict` and
`metadata_dict`, CCD from `lmdb/ccd/biomol_CCD.lmdb`) on each record and writes the `to_bytes` value. The schema is **A**:
a record has the shape of a `pdb/cif` value, with both assemblies of the CIF (`1_1_.`, `2_1_.`) in `assembly_dict`
(the example row says which assembly and which glycan unit a training sample takes).

**Our reader parses these CIFs, but the recipe returns an empty `assembly_dict` for them as they are.** The files are
written by AtomWorks and lack tables `pdb/cif` takes from the mirror's files. The adapter (in the script, to become a
reader `readers/reglyco.py`) derives them from what the CIF carries and invents nothing:

| Missing or different in the Re-Glyco CIF | What the adapter does |
|---|---|
| `_pdbx_poly_seq_scheme`, `_pdbx_nonpoly_scheme`, `_pdbx_branch_scheme`, `_pdbx_entity_nonpoly`, `_pdbx_entity_branch_list` | built from `_atom_site` (label / auth ids, residue names) and `_entity.type`; branched residues are numbered 1..n in chain order |
| `_entity_poly_seq` has no `hetero` | `n` (no heterogeneity is recorded) |
| `_chem_comp` lists the deposited entry's components only (no BMA, MAN, FUC, ...: 2,042 added over 700 records) | the missing ids are added by name; atoms and bonds come from the CCD |
| `_chem_comp_atom` / `_chem_comp_bond` lack `pdbx_stereo_config` | dropped, so the recipe takes the component from the CCD as it does for any CIF without them (filling the column would invent chemistry) |
| `_struct_conn` has `ptnr*_label_seq_id`, no `ptnr*_auth_seq_id`, and holds the sugar-sugar links of a branched chain | `auth` ids mapped from `_atom_site`; the intra-branch rows (6,337 in 700 records) move to `_pdbx_entity_branch_link`, as in a PDB mmCIF. The leaving atoms of those links are not in the CIF and are written `?` |
| no deposition date, no resolution, `_entry.id` is a temp file name | `_pdbx_database_status` and `_refine` are filled from `examples.parquet`; the id is the record key. `release_date` falls back to the deposition date (no revision history) |

**Check of the conversion** (all 700 records, 1,059 assemblies): for every assembly the number of atoms with finite
coordinates, and the sums of x, y, z, occupancy and B factor, equal those of the CIF's `_atom_site`. BioMol records also
carry every atom of a residue's CCD definition, the unobserved ones with NaN coordinates (7.1 M atoms out for 6.3 M
`_atom_site` rows), so the totals differ and only the finite atoms are compared. Imputed atoms stay (occupancy 0: 426 k
atoms), B factors stay (NaN where the CIF has `nan`). All polymer chains cover their `_entity_poly_seq` (the CIF imputes the
missing residues).

**What the record does not hold**, and where it goes:

| Dropped from the structure value | Why / where |
|---|---|
| per-atom AtomWorks columns of `_atom_site`: the PLIP `condition_*` / `mask_*` (hbond, hydrophobic, saltbridge, metal-coord, pistack, pication, cath-label), `nhyd`, `atom_array_index` | derived from the coordinates; no field for them in schema A. A model that conditions on PLIP contacts has to recompute them, or the stage keeps them in a sidecar |
| `_atomworks_imputation_stats` (QC scalars of the imputation) | not in `examples.parquet` either; the pass / fail flags and clash counts that are there stay in the sidecar |
| `_atomworks_cath` (CATH label per chain) | dropped; recoverable from the CIF in the raw download |
| `_link_restore` (a note on a deposited sugar link that was removed and rebuilt) | 313,306 records (10.5 %) carry one; the count `n_deposit_sugar_links_removed_with_sugar` stays in the sidecar |
| `_cell`, `_entry.author` / `date` / `time`, `_chem_comp` name / type | no use downstream |
| the 90 columns of `examples.parquet` (glycan, site, UniProt, sampling weights, QC flags, split) | **sidecar table** `metadata/reglyco_examples.parquet`, joined on `path` (key = file name without `.cif`); only `deposition_date` and `resolution` go into `metadata_dict` as well |

## Layouts on BioMol records

700 records, 5 entries, 108 groups. Shipped CIF: 159.2 MB. Converted (`to_bytes`): 119.7 MB (171 KB per record). A
delta is zstd with the serialised base record as a raw-content dictionary (window sized to base + record). The `to_bytes`
header names its arrays with random UUIDs; the codec renames them to `a0, a1, ...` before the delta (the decoded dict is
unchanged; without it A1 is 14-42 % larger). Every layout was decoded again and compared with the converted dict (arrays
and scalars equal): **700 of 700** in every row.

| Layout | Level | Size | Entry bases | Group bases | Conformer deltas | Encode / decode per record |
|---|---|---:|---:|---:|---:|---|
| as shipped (CIF, zstd) | - | 159.2 MB | | | | |
| one BioMol value per record | 6 | 119.7 MB | | | | |
| A1: entry base + a delta per record | 9 | 8.12 MB | 0.86 MB | | 10.4 KB each | 5.5 / 5.2 ms |
| A1 | 19 | 2.67 MB | 0.74 MB | | 2.8 KB each | 55 / 5.8 ms |
| A2: + group base (delta vs entry base) + conformer deltas vs the group base | 9 | 3.02 MB | 0.86 MB | 10.8 KB each | 1.8 KB each | 4.2 / 5.6 ms |
| **A2** | **19** | **1.99 MB** | 0.74 MB | **2.8 KB each** | **1.6 KB each** | 29 / 5.9 ms |

- **Level 19 matters for BioMol.** The group base (a different glycan: its atoms and edges change length, so every array
  after them shifts) is 10.8 KB at level 9 and 2.8 KB at level 19; the conformer delta (same glycan, new coordinates) is
  1.8 against 1.6 KB. For CIF text levels 6-12 were equal; here they are not. Level 12-15 were not measured.
- Decode 5-6 ms per record covers building the dictionary, the zstd decode and `np.load` of the arrays (about 170 per
  record); a loader that keeps the entry base's dictionary needs less. It is not the cost to worry about.
- Encode at level 19 is 29 ms per record (A2), 55 ms for A1.
- Not measured: B (one solid frame per group, which gave 26 MB on CIF text), other shards, entries with many glycans.

## Size for the whole set

Exact population from the Hub tables (`reglyco_extrapolate.py`); sizes extrapolated from the sample
(`reglyco_biomol_extrapolate.py`): the entry base is `r_base` times the shipped length of the entry's first record (BioMol:
0.61 at level 19, 0.71 at level 9), the delta means come from the sample, low / high are the smallest / largest per-entry
means of the five entries (indicative, not a confidence interval).

| | records | shipped CIF | BioMol, one value per record | **BioMol A2, level 19** | BioMol A2, level 9 | BioMol A1, level 19 |
|---|---:|---:|---:|---:|---:|---:|
| all | 2,992,448 | 1,573.6 GB | 1,183 GB (1,051-1,220) | **7.5 GB** (6.4-8.3) | 12.0 GB (9.2-14.1) | 10.3 GB (9.3-11.6) |
| training view | 2,307,945 | 1,192.1 GB | 897 GB (796-924) | **5.8 GB** (4.9-6.4) | 9.3 GB (7.1-10.9) | 8.0 GB (7.2-8.9) |

The training view is the one the configs use: deposited before 2024-01-01, `glycan_experiment_split == 'train'`, an
in-vocabulary glycan (2,307,945 rows, as the README says). Add 5-10 % for LMDB pages and keys.

**What the sample cannot tell.** Five entries are 0.04 % of the set, and their records are smaller than the population's
(1.4 MB against 3.6 MB on average). The group-base part is the one that could grow: it depends on how different two
glycans' atom counts are, which the five entries sample poorly. Run `LEVEL=19 PHASE=evaluate` of
`reglyco_biomol.sbatch` on a few shards of other size before fixing the schema.

## Build cost

- **Conversion**: 1.3 s of CPU per record (0.48 s reading the CIF, 0.81 s in the recipe), so **about 1,100 core-hours for
  all records (850 for the training view)**; on one 112-core node about 10 hours.
- **Delta encode at level 19**: 29 ms per record, about 25 core-hours.
- Both are small next to the 1.57 TB download; the conversion is the part to optimise if it matters: the protein context
  of an entry is the same in every record, so an entry's base could be converted once and a conformer produced by
  replacing the glycan's atoms. Not done.

## Fitting the delta codec to the DB contract

The loader should see only ordinary BioMol records. A delta value is not a `to_bytes` frame (a plain `from_bytes` fails on
it), so the DB needs its own codec and a reader; everything else of the contract stays.

- **Keys**: every record keeps its key (`<pdb>_<site>_<glycan>_conf<NNN>`, 2,992,448 keys). No reserved keys: **the bases
  are real records of the same DB**. The base of an entry is its first record by key, stored as a plain `to_bytes` value
  (so `from_bytes` reads it); the base of a group is its lowest conformer, stored as a delta against the entry base; every
  other conformer is a delta against its group base. Decoding takes at most two steps.
- **Value**: `to_bytes` frame for an entry base; otherwise `magic | ref_key | zstd frame (raw-content dictionary = the
  base's canonical serialised bytes)`. A new codec `biomol_zstd_delta` next to `biomol_zstd` in `build/schemas.py`, and a
  reader `read_record(env, key)` (and a cache of the entry's base) next to `from_bytes`; the loader calls that and gets the
  ordinary dict.
- **`.index.tsv` / `.meta.json`**: unchanged in meaning (key, stored size, counts); sizes are the delta sizes.
- **`release-check` and the schema validators**: the per-record check decodes a sample through `read_record` for DBs whose
  codec is `biomol_zstd_delta` and then runs the schema A validator; one new branch, keyed on the codec.
- **Sharding**: 2.99 M records are below the 5 M merge limit, so one LMDB of about 8 GB. A work item is one entry (all its
  records, since a delta needs the base), written with `explode_entries`-style output of already-compressed bytes.
- **Patching**: a record can only be patched together with the rest of its entry (a changed base invalidates its deltas).
  The patch unit is the entry.
- **Alternatives**: a sidecar DB for the bases would put two keys per record family and break the one-DB rule; a solid frame
  per group (B) has no per-record access and was 11 times worse on CIF text. A structure that stores the protein once and
  the glycans separately would change the schema and is not worth it at 7.5 GB.

## Recommendation

**BioMol schema A, delta layout A2 at zstd level 19, bases as real records of the same LMDB, a sidecar table for the
Hub's per-record metadata.** 7.5 GB (5.8 GB for the training view) instead of 1,574 GB, every record readable as an
ordinary BioMol dict through `read_record`. Build it after the download passed `verify --sha`: a `materialize` stage over
`blob/index.parquet` lists the entries and groups, a reader `readers/reglyco.py` holds the adapter above, a recipe
encodes the deltas and decodes each once against the converted value (the build fails on any mismatch).

Decision for the owner: whether to keep the 23 % of records outside the training view (about 1.7 GB; keeping them leaves
the holdout available), and whether a model needs the dropped PLIP masks and CATH labels (then a sidecar for them is part
of the stage).

## Result of the build

`lmdb/glycan/reglyco/reglyco.lmdb` was built on 2026-11-11/12 (`glycan/reglyco_worklist` -> `reglyco_entries` -> `reglyco`; 8 nodes, 5.5-7.4 h per node,
about 11,000 records per minute) from the sha256-verified download. Measured, not extrapolated:

| | estimated above | built |
|---|---:|---:|
| records | 2,992,448 | 2,992,448 (11,160 plain entry bases, 2,981,288 deltas) |
| size, all records | 7.5 GB (6.4-8.3) | **14.95 GB** (entry bases 2.18 GB, deltas 9.64 GB, the rest LMDB pages and keys) |
| shipped CIF for comparison | 1,574 GB | 1,574.6 GB, i.e. **105 times smaller** |
| conversion | 1.3 s per record | about 4 s per record (3,000-group sweep), roughly 3,300 core-hours |
| decode per record | 5-6 ms (sample) | cold (new reader, no cache) mean 51 ms, median 47, p95 99; warm (entry bases cached) mean 36 ms |

The estimate was low by a factor of two because the five sampled entries were smaller than the population: the records average 3.6 MB against
1.4 MB in the sample, and the entry bases and group-base deltas grow with the record. The build time was set by the same thing. Decode time is dominated by
rebuilding the arrays, not by zstd, so the entry-base cache saves less than the sample suggested (36 against 51 ms).

Checks of the finished DB (`scripts/analysis/reglyco_check.py verify`): every one of the 2,992,448 records decodes to a BioMol dict with a non-empty
`assembly_dict` and a deposition date (0 failures, 31 min on 60 workers); 5,000 random records (4,974 of them delta values), each read with a cold reader, equal a fresh
conversion of their CIF (canonical serialised bytes; 0 mismatches). A sweep of one record in each of 3,000 random (entry, site, glycan) groups converted without
failure after the adapter learned lone sugars (below). The pipeline's index step wrote `codec: biomol_zstd` into the DB's meta.json; it was set to
`biomol_zstd_delta` by hand, and `release-check` now fails a delta-coded DB whose meta.json does not say so. Since 2026-11-13 the build passes the schema's codec to the pipeline (`--codec`), so a rebuild writes it itself (the DataCooker change is a local commit until that repository is pushed).

Two things that went wrong, for the next build: (1) a one-residue glycan is typed `branched` by AtomWorks without a sugar-sugar link, which the `pdb/cif` recipe cannot read (the entry
failed); the adapter now types such an entity `non-polymer`, as a PDB mmCIF does. (2) For an LMDB rebuild the engine runs the first 10 keys of each shard serially before
the parallel run, and here those are the ten largest entries (about 20,000 records, a day of one core): `test_run: false` in the config. The engine's admission gate also starts at one task and
grows by one per finished wave, which left most cores idle for hours; `DATACOOKER_ADMISSION_START=32` (a local DataCooker commit, see
[development](development.md)) starts at 32 per node. The memory gate still applies; the run peaked near 110 GB of 503 GB per node.

## Layouts on raw CIF text (the earlier measurement)

Before the BioMol requirement was applied, the same shard was measured on the CIF text itself. Kept as the lower bound
of what the information is worth; it is **not** the recommended form, because it stores CIF bytes, not BioMol records.

700 records, 5 entries, 108 groups. Shipped: 159.23 MB (227 KB per record); decompressed 972.6 MB. Every layout was
decoded again and compared with the original bytes and the index's `cif_sha256`: all 700 records match.

| Layout | What is stored | Size | Per record (A1/A2: the delta part) | Decode |
|---|---|---:|---:|---|
| shipped | one zstd frame per record | 159.2 MB | 227 KB | - |
| per-record, level 6 | the same, recompressed | 164.1 MB | 234 KB | - |
| B, level 12 | one solid zstd frame per group | 26.1 MB | 37 KB | 2 ms per group |
| A1, level 9 | entry base + a delta per record | 2.84 MB | 2.6 KB | 0.1 ms |
| A2, level 9 | entry base + group base + conformer deltas | 2.31 MB | 1.7 KB per conformer | 0.1 ms |
| A2, level 19 | as A2 | 1.99 MB | 1.3 KB per conformer | 0.1 ms |
| A3 (sample) | entry base + explicit line diff, zstd 19 | - | 2.1 KB | 1 ms |

Extrapolated: A2 level 9 **8.4 GB** (7.2-9.4; training view 6.5 GB), level 19 7.0 GB (6.2-7.8; 5.5 GB). The delta part of
a CIF text record is about 2 KB, the information floor of the 80 changed glycan lines. Level 3 does not work (the match
finder does not find the 1.9 MB base). A "delta" there is zstd with the base record as a raw-content dictionary; build
cost on this form is about ten core-hours, because nothing is parsed. Code: `reglyco_compact.py`,
`reglyco_extrapolate.py`.

## Reproduce

```bash
source scripts/cluster/cssb.env
sbatch scripts/analysis/reglyco_biomol.sbatch                      # convert rs_00005 + A1/A2 at level 9 (2 min)
LEVEL=19 PHASE=evaluate sbatch scripts/analysis/reglyco_biomol.sbatch
sbatch scripts/analysis/reglyco_biomol_extrapolate.sbatch          # population + sizes -> extrapolation.json (10 s)
sbatch scripts/analysis/reglyco_compact.sbatch                     # the raw-CIF layouts (3 min)
sbatch scripts/analysis/reglyco_extrapolate.sbatch                 # raw-CIF population + sizes
```

`examples.parquet` was fetched separately to `materials/intermediate/reglyco_compact/hub/`; the scripts read it from
there.
