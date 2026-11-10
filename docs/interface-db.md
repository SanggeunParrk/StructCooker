# Interface-unit DB (PDB)

A DB of **interface units**: groups of 2..N molecules of one PDB bioassembly that are in contact,
stored with all their heavy-atom coordinates **and the surrounding waters**. A unit is the sample
for models that need "a binder and its target" (antibody + antigen, protein + dsDNA, protein +
ligand + ions, ...). Voxelisation or cropping is left to the consumer.

Code: `src/structcooker/instructions/interface/` · thresholds: `src/structcooker/instructions/interface/params/` ·
stages: `db_configs/interface/` (set manifest `db_configs/MANIFEST_interface.yaml`, part of `MANIFEST_all.yaml`) ·
recipes: `src/structcooker/workflows/` · developer tools: `scripts/interface/` · tests: `tests/interface/`. The DB is a normal
part of the release DAG: `structcooker build interface/<stage>`, `build-all`, and `e2e-sample` build it like any other.

## What is a node, what is a unit

A **node** is one of

| kind | what | key |
|---|---|---|
| `protein` | a polypeptide chain | `P:<sequence>` (canonical parent letters) |
| `na` | a DNA / RNA / hybrid chain | `N:<sequence>` |
| `ligand` | a CCD molecule, or a glycan chain; a ligand covalently bonded (`struct_conn` `covale`) to a polymer is merged into that polymer's node | `L:<ccd>` or `L:<ccd-ccd>` |
| `ion` | a single-atom non-polymer (Mg, Zn, Ca, ...) | `I:<ccd>` |

Water is never a node.

**Edges** (all thresholds in `src/structcooker/instructions/interface/params/interface.yaml`):

- polymer-polymer: heavy atoms < 5.0 A apart, and >= 5 contacting residues on *both* sides
- protein-ligand: >= 3 ligand heavy atoms within 5.0 A and >= 2 contacting polymer residues
- ion: N / O / S donors of a polymer within 3.0 A (MetalPDB and InterMetalDB use the same distance); an
  ion starts a unit only when it bridges >= 2 nodes (an ion on one chain is not an interface)
- NA-NA and NA-ligand are not edges (extension point: `unit.allowed_pairs`)

**Units.** Every significant edge seeds one unit. R = the contacting residues of both sides of the seed
edge; every other node that has a significant edge touching R joins (one hop only, so a complex does not
percolate): >= 3 of its residues touch R for a polymer, one contact for an ion or ligand. Caps: 6 polymers,
4 ligands, 8 ions per unit; a unit that hit a cap is flagged `truncated`. Units with the same node set are merged
(`n_seed_edges`), units equal in node sequences and contact size are deduplicated within an entry (symmetry copies
and repeated assemblies), assemblies with more than 60 polymer chains are skipped.

**Forced partners** (bypass the caps, up to 10 polymers):

- *antibody / TCR*: chains are typed from sequence with ANARCI (IMGT): H, K, L, TCR A, B, G, D. A unit that
  holds one gets its best-contacting complementary chain (H with K/L, A with B, G with D). A bare H-L unit
  is dropped when a larger unit contains the pair. Nanobodies and scFvs are single chains.
- *double-stranded NA*: two NA chains are a duplex when they share >= 4 base pairs (C1'-C1' 8-12 A and >= 2
  base N/O pairs within 3.5 A; mismatches, wobble, Hoogsteen count; DNA, RNA and hybrids alike). A unit that holds
  a strand gets every strand it base-pairs with, whether or not the partner touches the protein.

**Ligand exclusion** (`src/structcooker/instructions/interface/params/ligand_exclusion.yaml`, categories switched in the config): the AF3
crystallisation aids (SO4, GOL, EDO, PO4, ACT, PEG, DMS, TRS, PGE, PG4, FMT, EPE, MPD, MES, CD, IOD; transcribed
from the AF3 SI from memory, not yet checked against the PDF), further common additives, water, monovalent
ions / halides, detergents and lipids. Other ions are nodes. Ligands below 6 heavy atoms are dropped.

**Water.** Stored as context: every water O within 4.0 A of any unit atom. `bridging` = within 3.5 A of atoms of
two or more different unit nodes. Waters never enter edges or residue sets. Water is absent from `cif_pdb_attached`
(it has no sequence, so it cannot get a seq_id and would break the chain set the FASTA / cluster tables are built
on); the interface DB reads the base `cif_pdb` instead, which keeps it.

**Other rules.** First model of each assembly (NMR entries: model 1); altloc: the no-altloc `.` entry, else the first;
modified residues map to their parent through the CCD canonical one-letter code (unknown -> `X` / `N`); atoms with no
coordinates are not stored; every bioassembly of an entry is processed, and units repeated across assemblies are kept once.

## Stages

Five stages in `db_configs/interface/`, in this order (inputs and outputs in [stages](stages.md#interface-units)):

| Stage | Op | What it does | Code |
|---|---|---|---|
| `interface/antibody_cache` | materialize | ANARCI (IMGT) chain type of every unique PDB protein sequence of length 70-600 -> `antibody_cache.tsv` | `antibody.annotate_antibodies` |
| `interface/interface_pdb` | rebuild, `explode_entries`, schema `J` | every `cif_pdb` entry -> its unit records | `units.explode_units` (reader adapter), `units.build_unit_record` (recipe step), `record.unit_to_record` |
| `interface/manifest` | extract | one row per unit read back from the DB -> `manifest.tsv` | `record.manifest_row`, `tables.write_manifest` |
| `interface/manifest_full` | materialize | manifest + each protein node's seq_id and 40% cluster, `cluster_key40`, `including_NA_entry` -> `manifest_full.tsv`; statistics -> `report.md` | `tables.write_manifest_full` |
| `interface/split` | materialize | date split with leakage flags -> `split.tsv`, `split_summary.md` | `tables.write_split` |

`interface_pdb` is the 1:N rebuild (like `chain/cif_chain`): the engine writes one record per sub-entry the reader adapter
returns. The adapter has to see the whole entry (first model of each assembly, units repeated across assemblies dropped,
assemblies with more than 60 polymer chains skipped), so it finds the units and returns one sub-entry per unit, keyed
`<assembly>_<model>_<index>`; the recipe builds the record, the engine keys it `<pdb id>_<assembly>_<model>_<index>`. The
adapter reads the ANARCI cache through `$INTERFACE_ANTIBODY_CACHE`, which the stage config declares in
`runtime_environment` (so the cache is a fingerprinted, pinned input of the stage); a sequence missing from the cache is
classified live. Protein nodes are matched to sequences by `md5("P:<sequence>")` over the PDB protein fasta.
The thresholds sit next to the code, not in `db_configs/`: a release runs a frozen copy of `src/`, and changing them
changes the code the release refuses to resume with.

## Where it lives

```
$BIOMOL_ROOT/lmdb/interface/pdb/
  interface_pdb.lmdb             the units: 850,925 records (a merged LMDB below 5M records, else a shard collection)
  interface_pdb.lmdb.index.tsv / .meta.json   written by the engine; what release-docs lists
  manifest.tsv / manifest_full.tsv   one row per unit; _full adds seq_ids and 40% clusters
  split.tsv, split_summary.md    date split with leakage flags
  report.md                      counts and distributions
$BIOMOL_ROOT/materials/intermediate/interface/antibody_cache.tsv   ANARCI type of every PDB protein sequence
```

The release published on 2026-10-28 was built by the earlier driver and is a 141-shard collection with a `shard` column in
the manifests and a per-entry `done.tsv`; the stages produce the same records (checked) but one merged LMDB, no `shard`
column and no `done.tsv` (entries without a unit are counted from `cif_metadata` in `report.md`).

Key = `<pdb>_<assembly>_<model>_<index>` (e.g. `1tsr_1_1_0`). Values use the repository codec
(`codecs.to_bytes`: zstd + numpy arrays). Records are self-contained: no lookup into `cif_pdb` is needed.

```python
from structcooker.instructions.interface.db import InterfaceDB
db = InterfaceDB("/data/shared/cssb_data/BioMol_clean/lmdb/interface/pdb")
rec = db.get("1tsr_1_1_0")        # reads a merged LMDB, a shard collection or bare shards/
rec["meta"], rec["nodes"], rec["atoms"], rec["residues"], rec["edges"], rec["duplex"], rec["water"]
```

Schema (v1; the docstring of `record.py` is authoritative):

- `meta`: pdb_id, assembly_id, model_id, alt_id, unit_id, resolution, release_date, deposition_date, signature, truncated, n_seed_edges
- `nodes` (per node): kind, key, hash, chain_id (asym_oper), auth_asym_id, entity_id, ab_class
- `atoms`: xyz (original PDB frame), element, atom_name, node, res, b_factor, occupancy, contact
- `residues`: chem_comp_id, parent, cif_idx (position in the node's sequence), auth_idx (author number, may be `95.A`), node, interface
- `edges`: a, b, n_atom_pairs, n_res_a, n_res_b, min_dist, center, res_pair_ptr / res_pairs
- `duplex`: a, b, n_bp
- `water`: xyz, chain_id, auth_idx, bridging

Why self-contained: a unit holds 1 to 10 chains of an assembly that may have 60, so referencing the source record
would mean decoding up to a 1 GB record per sample; zstd takes most of the overlap between a complex's units out and
the consumer wants random access.

Not stored: buried SASA, the experimental method (X-ray / EM / NMR is not in `cif_pdb`).

## Split

Date based, as in AF3, and per PDB entry (every unit of an entry lands in the same split). `split.tsv`:

| split | rule | units | entries |
|---|---|---:|---:|
| train | released before 2021-09-30, resolution <= 9 A (the rule of `train_20210930`) | 488,379 | 134,324 |
| valid | 2021-09-30 .. 2023-01-12 (AF3's PDB evaluation window), resolution <= 4.5 A | 86,897 | 14,048 |
| test | after 2023-01-12, resolution <= 4.5 A | 248,054 | 40,566 |
| none | no release date, no resolution (NMR), or above the cut | 27,595 | |

valid and test units carry leakage flags against the training set (nothing is dropped): `n_prot_c40_in_train`,
`novel_protein` (no protein of the unit has its 40% cluster in train), `novel_pair` (the unit's set of clusters
does not occur in train), `ligands_seen_in_train`. Of the valid units, 8,184 have no protein cluster seen in train and
22,179 have a new cluster combination; filter on these for a strict evaluation set.

## Reproduce

```bash
structcooker build interface/antibody_cache    # once: ANARCI over every unique protein sequence
structcooker build interface/interface_pdb     # the units: rebuild of cif_pdb, sized by the engine from the entry sizes
structcooker build interface/manifest          # manifest.tsv
structcooker build interface/manifest_full     # manifest_full.tsv, report.md
structcooker build interface/split             # split.tsv, split_summary.md
# or all of it with its upstream: structcooker build-all --manifest db_configs/MANIFEST_interface.yaml --workdir <run dir>
pixi run python scripts/interface/export_unit.py $BIOMOL_ROOT/lmdb/interface/pdb /path/to/cifs 1tsr_1_1_0
pixi run python scripts/interface/validate_db.py $BIOMOL_ROOT/lmdb/interface/pdb 1000   # recompute contacts / flags / waters (a check, not a stage)
```

`structcooker build NAME --show` prints the resolved invocation, `--dry-run` writes the sbatch without submitting. A
change of the thresholds (`params/interface.yaml`) is a code change: rebuild `interface/interface_pdb` and everything after
it into a new tree. `structcooker e2e-sample` builds these stages with the rest of the DAG on a sample and requires its unit
records to equal the release's, key by key. Peak memory is ~16x the compressed record (5j7v, a 1.3 GB entry, needed 21 GB
in the earlier driver); the engine plans memory from the entry size (schema `J`, E = 100). `scripts/interface/run_stats.py`
and `run_sample.py` are the sampling tools used to choose the thresholds; `make_viewer.py` builds a 3D viewer page of
example units (`pick_examples.py` selects them).

`export_unit.py` writes mmCIF with B_iso as a marker: 100 interface residue, 50 other contact atom, 20 water (80
bridging), 0 the rest (`--keep-bfactor` keeps the real ones).

## Checks

`tests/interface/` (synthetic geometry: contact search against brute force, joining, caps, forced partners,
duplex, ions, waters; the unit adapter / recipe / explode chain, the antibody cache, manifest, report and split). `validate_db.py` recomputes every edge's min distance, atom-pair count (within 1e-3 A of the
cutoff: coordinates are stored as float32), the contact flags, water distances and the bridging flag from the stored
coordinates alone.
