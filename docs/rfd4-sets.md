# RFD4-Proteina structure sets as BioMol records

Seven components of the Hugging Face dataset `nvidia/RFD4-Proteina` (Hub commit `aa86df1e`, downloaded raw: [raw-materials.md](raw-materials.md#rfd4_proteina)) are converted to
ordinary BioMol records, one DB per component under `lmdb/rfd4/<set>/<set>.lmdb`. The eighth downloaded component, Re-Glyco, has its own
stages ([reglyco-compact.md](reglyco-compact.md)). The PDB stays one DB; these are the non-PDB sets a data loader mixes in.

| set (`rfd4/<set>`) | Hub component | records | what a record is | key |
|---|---|---:|---|---|
| `omol` | `omol` | 53,664,160 | an OMol molecule or complex (DFT geometry; small molecules, metal complexes, peptide fragments) | zero-padded `lmdb_idx` (`000000007`) |
| `rna_monomer` | `rna_monomer` | 371,017 | one predicted RNA chain (Protenix-v2) | `example_id` (`URS000080E21C_32630`) |
| `na_complex` | `na_complex` | 20,261 | a predicted protein–DNA complex (ESMFold2) | `example_id` |
| `antibody_distill` | `antibody_interfaces_distill` | 75,576 | a predicted antibody–antigen complex (ESMFold2) | `example_id` (`aae_000022`) |
| `tcr_distill` | `tcr_interfaces_distill` | 3,941 | a predicted TCR / peptide–MHC complex (AF2) | file stem of `path` (`tcr_10`) |
| `synthetic_ncaa` | `synthetic_ncaa` | 1,463,222 | one peptide conformer with PTMs / NCAAs (ESMFold2) | `example_id` |
| `ncaa_scaffolds` | `ncaa_scaffolds` | 288,390 | one NCAA scaffold peptide conformer (ESMFold2) | `example_id` |

The record keys, the sha256 column and the default date of each set are in one table, `COMPONENTS` in
`src/structcooker/instructions/transforms/rfd4_worklist.py`.

## One format

A record is the value of the PDB's own build, schema A's shape `{assembly_dict, metadata_dict}` (registered as schema **M**: the same, with a free record id as the key).
It is made by the **same recipe** as `pdb/cif` (`workflows/ingest/cif.py`); nothing downstream needs to know which set a record came from. No special codec.
`metadata_dict` has `id` (the key), `deposition_date` (see below) and `release_date`; `resolution` is `None` (predicted or simulated structures).
Everything else the Hub says about a record (scores, clusters, sampling flags, the per-set interface annotations) is **not** in the value: it is in the sidecar
`metadata/rfd4_<set>_examples.parquet` (the Hub's `examples.parquet` plus a `key` column; join on `key`). The TCR set has two rows per record (one per
receptor / pMHC unit pair); the others have one.

## Stages

```
rfd4/<set>_worklist  (materialize)  blob/index.parquet + examples.parquet  ->  materials/intermediate/rfd4/<set>/worklist.lmdb  and the sidecar
rfd4/<set>           (LMDB rebuild, 1:1)  worklist  ->  lmdb/rfd4/<set>/<set>.lmdb    needs ccd/ccd
```

The worklist is one tiny row per record: `{key, shard, offset, length, sha256, date[, sequence]}`; `shard` is relative to the component's `blob/`
(the shards of the 2026-10 reannotation sit beside `blob/data/`, so the stage finds each by looking, once per distinct shard). The rebuild's reader
adapter (`readers.rfd4.rfd4_record_to_cif_inputs`) reads the record's frame, checks the decompressed CIF's sha256 against the index, and hands the recipe its inputs.
A record that fails (bad frame, a recipe error) is counted and skipped, as in `pdb/cif`; `scripts/rfd4/validate_db.py` lists the worklist keys a DB lacks.

```bash
source scripts/cluster/cssb.env
pixi run structcooker build rfd4/<set>_worklist
pixi run structcooker build rfd4/<set>          # DATACOOKER_ADMISSION_START=32 is worth setting for omol
```

## What a thin CIF needed (`readers/rfd4.py`)

The Hub's CIFs are written by AtomWorks from prediction and simulation pipelines and are thinner than a PDB mmCIF: some have only `_atom_site` and a partial
`_chem_comp_bond`; four of seven have no `_entity` at all. The adapter derives what the recipe reads, **only from what `_atom_site` states**:

* `_entity`: the ids are the atom_site's `label_entity_id`; the type is `polymer` for a chain holding a standard residue (or several residues numbered in
  sequence), `water` for HOH, else `non-polymer`. `_entity_poly` (type from the residue classes: protein, DNA, RNA), `_entity_poly_seq` and `_chem_comp` likewise.
  A table the record has is kept; a column of `_entity_poly` it lacks (`pdbx_seq_one_letter_code*`) is filled from its own sequence table.
* the canonical sequence letter of a residue is its standard letter; a residue that is not standard is `X` (`N` in a nucleic acid), as for any unknown
  component. Where the dataset states a canonical sequence for a record with one polymer chain (`synthetic_ncaa`, `ncaa_scaffolds`: column `sequence`)
  and its length equals the chain's, that sequence is used, so a PTM keeps its parent letter.
* the scheme tables, the branched-glycan tables and the `_struct_conn` fields come from `readers.reglyco.adapt_annotated_cif` (the same adapter Re-Glyco uses).
* `occupancy` (1.0) and `B_iso_or_equiv` (0.0) are filled when a prediction has none; `deposition_date` is the record's own `_entry.date` when it is a date,
  else the component's release date (`COMPONENTS[...].default_date`). That is a release date, not a deposition: no cutoff by date means anything inside these sets.
* components: a CIF's `_chem_comp_bond` without atoms is dropped (the recipe takes the component from the CCD, as for any CIF without `_chem_comp_atom`).
  A record with `_chem_comp_atom` lists every component it uses, standard ones included, without the CCD's stereo flags: the CCD components are dropped from the
  record's table and come from the CCD whole; what stays are the record's **own components** (an OMol molecule is named `A:1`, `B:1`, ... and a peptide
  cap `A:N` / `A:C`), which no CCD defines.

`transforms.cif.compare_chem_comp` (the recipe's component step) therefore learned one thing: a component the CCD does not know but the record defines is
the record's own definition (before, the lookup raised `KeyError` and the entry was lost). Its atom table gets `model_xyz` = `?` (a CCD component carries ideal
coordinates, and the atom tables of one record must have the same features). This also applies to a PDB entry with a component newer than the CCD.

## Not carried

* Per-residue / per-atom annotations of the CIFs that are not BioMol features: `ig_region`, `cdr` (antibody and TCR sets), the PLIP `condition_*` / `mask_*`
  columns, CATH labels, repair masks, `atom_array_index`. The BioMol record has coordinates, bonds and entities. The per-record summaries
  (`cdr_annotations`, `ig_chains`, `interface_cluster`, ...) are in the sidecars; the per-residue `ig_region` / `cdr` are in `metadata/rfd4_antibody_distill_residues.parquet` and `metadata/rfd4_tcr_distill_residues.parquet` (one row per residue: key, label / auth chain and residue ids, component, `ig_region`[, `cdr`]; `scripts/rfd4/residue_annotations.py`).
* OMol hydrogens: its CIFs carry them (about half of an OMol molecule's atoms; 144 atoms in the CIF of `009018428`, 73 in its record), the PDB recipe drops hydrogens from every
  record, so an OMol record has none, and the Hub's `n_atoms` (which counts them) is not comparable with the record's atom count. A residue's missing standard atoms are
  listed from the CCD, as in any PDB record. The CIFs in the shards keep the hydrogens if a consumer needs them.
* OMol: its `examples.parquet` has 53,485,498 rows for 53,664,160 index records (the Hub removed 836 failed records from the index *and* sampling rows; the rest of the
  difference is index records without a sampling row); both are in the DB / the sidecar as they are. The Hub's own warning stands: OMol's source metadata
  identifiers, categories, charges and energies are not reconciled with the CIF block `lmdb_idx` retrieves.

## Checks

`scripts/rfd4/probe.py` (sbatch) converts a random sample of every set and prints failures and ms per record; `scripts/rfd4/validate_db.py` (sbatch) decodes every
record (or a sample) of a built DB, checks it has an assembly and a date, compares its residue counts with the sidecar (`n_standard_amino_acid_residues`,
`n_rna_residues`, `n_dna_residues`; the atom count may not be smaller than the Hub's: a BioMol residue lists the CCD's atoms, resolved or not) and lists the
worklist keys the DB lacks. Tests: `tests/rfd4/`, `tests/instructions/transforms/test_ccd_atom_order.py`.

## Result

Built 2026-11-12 from Hub commit `aa86df1e`, one `structcooker build` per stage (`DATACOOKER_ADMISSION_START=32` for omol), 8 nodes each.

| set | records in the DB / on the Hub | size on disk | converted by |
|---|---:|---:|---|
| `tcr_distill` | 3,941 / 3,941 | 232 MB | seconds |
| `na_complex` | 20,261 / 20,261 | 853 MB | seconds |
| `antibody_distill` | 75,576 / 75,576 | 9.2 GB | a few minutes |
| `rna_monomer` | 371,017 / 371,017 | 28 GB | a few minutes |
| `synthetic_ncaa` | 1,463,222 / 1,463,222 | 16 GB | a few minutes |
| `ncaa_scaffolds` | 288,390 / 288,390 | 8.8 GB | a few minutes |
| `omol` | 53,656,064 / 53,664,160 | 216 GB (8 shards, `omol.lmdb/shards/`) | about 1 h |

**OMol: 8,096 records are not in the DB** (0.015 %), listed with the reason in `metadata/rfd4_omol_unconvertible.tsv`: 8,042 name a molecule (`A:1`, `B:1`, a peptide cap
`C:N`, ...) that is in no CCD and has no `_chem_comp_atom` in its own CIF, so its elements and bonds are not stated anywhere (they fail with that sentence;
inferring bonds from the coordinates would invent chemistry), and 54 have an empty CIF. `patch` of the same keys with the final code fails the same 8,096.

Checks (all 2026-11-12):

* every worklist key but those 8,096 is in its DB (`validate_db.py`: *lacking 0* for the six others), and a decoded sample of each (40,000; antibody 8,000; tcr all; omol 50,000) has an assembly
  and a date, with residue counts equal to the sidecar's: protein / RNA / DNA residue counts match in every sampled record once the Hub's "non-standard" and "other" residues are
  added to the protein count (the NCAA sets count a modified residue as *other*); the atom count of a record is never below the Hub's (it is above it: the residue lists the CCD's atoms). OMol is the exception on both: its 1,890 of 50,000 sampled records differ in protein residues (peptide fragments with custom cap residues, `A:N` / `A:C`, that the Hub does not count as standard amino acids), and its atom counts are not comparable (hydrogens, above);
* the sampled end-to-end run (`e2e-sample prepare ... --count rfd4=40`, `run` and `verify --manifest db_configs/MANIFEST_rfd4.yaml`; 15 stages complete): 40 of 40 sampled records of each of the
  seven sets are identical to the release's, and every sampled sidecar row's key is in the release's sidecar;
* `release-check` (127 stages): every `rfd4/*` stage passes (entry counts agree with `.meta.json` and `.index.tsv`, 200 sampled records decode); the 10 stages it reports as problems are `design/*`, `metadata/ab_seq` and `metadata/cdr_region`, which are not part of this.

Two things the sampled run caught, now fixed: an output key of a materialize stage must start with `out_` (otherwise the DAG sees it as an input that changed), and a stage cannot declare a directory
as a pinned reader input (`runtime_environment` takes files), so the shards are found from `BIOMOL_ROOT` and the worklist row's `set`. The second also applied to `glycan/reglyco`, whose
config no longer declares them (`REGLYCO_BLOB_DIR` still overrides).
