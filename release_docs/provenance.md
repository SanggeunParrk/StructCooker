# Provenance

How this tree came to be, and what differs from a from-scratch build. Per-database counts
are in [databases.md](databases.md); the full build history is in StructCooker
(`db/STATUS.md`, `docs/history/`).

## Sources

| source | where | version |
|---|---|---|
| PDB mmCIF | `materials/raw/cif` | wwPDB rsync mirror (`structcooker download mmcif`) |
| CCD | `materials/raw/ccd` | wwPDB `components.cif`; split into `materials/intermediate/ccd/components` |
| OpenFold3 distillation | `materials/raw/openfold_distillation` | `s3://openfold3-data`: short/long monomers (released relaxed models `best_structure_relaxed.pdb`, preprocessed MSAs/templates), RNA set, disordered set |
| Teddymer | `materials/raw/teddymer` | non-singleton representatives passing the source paper's interface filter (510,454) |
| AFDB multimer | `materials/raw/afdb_multimer` | EBI `alphafold/collaborations/nvda`: high-confidence homodimers, heterodimers, release MSAs |
| SabDab | `materials/raw/sabdab` | antibody summary used to split antibody chains when clustering |

## Supplied inputs (`materials/reference`)

These cannot be derived from the downloads; each pins part of the result to the release.

| input | fixes |
|---|---|
| `seq_id_map_seed.tsv` | the seq_id numbering (this release's map); new sequences are appended after it |
| `msa/` | the PDB protein/RNA MSAs (HHblits, searched once) |
| `msa_teddymer/` | the Teddymer MSAs (MMseqs2 vs UniRef30 2302, searched once) |
| `openfold_distillation/{fasta,disordered_set}` | the OFD sequence set and the disordered work lists |
| `signalp/` | signal-peptide predictions used to trim training structures (7-digit legacy ids, matched by counter) |
| `reference_ccd.lmdb` | 24 historical CCD definitions older entries need |
| `reference_chain_selection.lmdb` | the per-chain model choice on occupancy ties |

## Known gaps and differences

- **PDB templates** (`lmdb/chain/template/*`) are the 2026-07-30 build. The chain DB gained
  293 entries (51,895 chains) on 2026-09-30, which are therefore not template candidates.
  Templates are restricted to hits released by 2021-09-30.
- **Failed inputs**: 5 PDB entries the parser cannot map and 840 OpenFold `alignment.npz`
  files that are empty or corrupt in the release (long 667, short 173) are absent from
  their databases (StructCooker `db/release-exceptions.json`).
- **OFD long**: 80 entries have no model in the OpenFold release; 14 entries whose
  `structure.npz` is truncated use recovered copies (`materials/intermediate/openfold_distillation/long_recovered`).
- **AFDB heterodimers**: 6 mmCIFs whose `_entity_poly` contradicts `_entity_poly_seq` are
  read from repaired copies (`materials/intermediate/afdb_multimer/heterodimer/repaired`);
  one homodimer whose MSA the release lacks uses an MMseqs2 MSA we generated
  (`materials/intermediate/afdb_multimer/msa_generated`).
- **pLDDT**: OpenFold publishes none for the RNA and disordered sets; their b-factors are 0.
- **MSA/template counts** are per distinct sequence, not per structure.

## Checks

- `structcooker release-check`: every database present at its path, self-contained,
  counts consistent with its sidecars, a sample of records decoding ([check.md](check.md)).
- `structcooker e2e-sample`: the whole build DAG re-run on a sample of every source;
  structures, MSAs and the CCD rebuild record-for-record identical to this tree.
