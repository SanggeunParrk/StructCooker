# Interface-unit DB (PDB)

A DB of **interface units**: groups of 2..N molecules of one PDB bioassembly that are in contact,
stored with all their heavy-atom coordinates **and the surrounding waters**. A unit is the sample
for models that need "a binder and its target" (antibody + antigen, protein + dsDNA, protein +
ligand + ions, ...). Voxelisation or cropping is left to the consumer.

Code: `src/structcooker/interface/` · thresholds: `configs/interface/` · scripts: `scripts/interface/` ·
tests: `tests/interface/` · commands: `structcooker interface`. Outside the `build-all` DAG: one input record fans out
into many unit records, so it has its own driver.

## What is a node, what is a unit

A **node** is one of

| kind | what | key |
|---|---|---|
| `protein` | a polypeptide chain | `P:<sequence>` (canonical parent letters) |
| `na` | a DNA / RNA / hybrid chain | `N:<sequence>` |
| `ligand` | a CCD molecule, or a glycan chain; a ligand covalently bonded (`struct_conn` `covale`) to a polymer is merged into that polymer's node | `L:<ccd>` or `L:<ccd-ccd>` |
| `ion` | a single-atom non-polymer (Mg, Zn, Ca, ...) | `I:<ccd>` |

Water is never a node.

**Edges** (all thresholds in `configs/interface/interface.yaml`):

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

**Ligand exclusion** (`configs/interface/ligand_exclusion.yaml`, categories switched in the config): the AF3
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

## Where it lives

```
$BIOMOL_ROOT/lmdb/interface/pdb/
  interface_pdb.lmdb/            shard collection (shards.json + shards/shard_NNN.lmdb), 850,925 records, 153 GB
  interface_pdb.lmdb.meta.json   what release-docs lists
  manifest.tsv / manifest_full.tsv   one row per unit; _full adds seq_ids and 30% / 40% clusters
  split.tsv, split_summary.md    date split with leakage flags
  report.md                      counts and distributions
  done.tsv, failed.tsv           per-entry build log (failed.tsv is empty)
$BIOMOL_ROOT/materials/intermediate/interface/antibody_cache.tsv   ANARCI type of every PDB protein sequence
```

Key = `<pdb>_<assembly>_<model>_<index>` (e.g. `1tsr_1_1_0`). Values use the repository codec
(`codecs.to_bytes`: zstd + numpy arrays). Records are self-contained: no lookup into `cif_pdb` is needed.

```python
from structcooker.interface.db import InterfaceDB
db = InterfaceDB("/data/shared/cssb_data/BioMol_clean/lmdb/interface/pdb")
rec = db.get("1tsr_1_1_0")        # pass shard=<manifest shard> to skip the search over shards
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

valid and test units carry leakage flags against the training set (nothing is dropped): `n_prot_c30_in_train`,
`novel_protein` (no protein of the unit has its 30% cluster in train), `novel_pair` (the unit's set of clusters
does not occur in train), `ligands_seen_in_train`. Of the valid units, 6,590 have no protein cluster seen in train and
20,360 have a new cluster combination; filter on these for a strict evaluation set.

## Reproduce

```bash
pixi run structcooker interface annotate     # once: ANARCI over every unique protein sequence
pixi run structcooker interface build        # whole PDB, 3 size tiers on one node (~30 min); rerun to resume
pixi run structcooker interface finalize     # manifest_full.tsv, report.md, publish the shard collection + meta.json
pixi run structcooker interface split        # split.tsv
pixi run structcooker interface validate --n 1000     # recompute contacts / flags / waters from the stored coordinates
pixi run python scripts/interface/export_unit.py $BIOMOL_ROOT/lmdb/interface/pdb /path/to/cifs 1tsr_1_1_0
```

The commands only sbatch (`--dry-run` prints the command); logs go to `logs/interface/`. Thresholds live in
`configs/interface/`; change them and rebuild into a new directory (a finalized directory refuses `build`). Before a
full build try `scripts/interface/build_db.py --limit 300 --out <scratch>`. Peak memory is ~16x the compressed record
(5j7v, 1.3 GB, needs 21 GB), which is why the tiers use 104 / 24 / 3 workers. `scripts/interface/run_stats.py` and
`run_sample.py` are the sampling tools used to choose the thresholds; `make_viewer.py` builds a 3D viewer page of
example units (`pick_examples.py` selects them).

`export_unit.py` writes mmCIF with B_iso as a marker: 100 interface residue, 50 other contact atom, 20 water (80
bridging), 0 the rest (`--keep-bfactor` keeps the real ones).

## Checks

`tests/interface/` (synthetic geometry: contact search against brute force, joining, caps, forced partners,
duplex, ions, waters). `validate_db.py` recomputes every edge's min distance, atom-pair count (within 1e-3 A of the
cutoff: coordinates are stored as float32), the contact flags, water distances and the bridging flag from the stored
coordinates alone.
