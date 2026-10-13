# Exports

Export workflows pull artifacts out of the LMDBs. Most read an LMDB through a
recipe (`extract` op — the config has `db_path` + `extract_recipe_path`); a few
derive tables from existing files (`materialize` op). `structcooker build`
infers the op and runs it through `datacooker.cli.workflow` as one SLURM job.

## FASTA

Extract all chain sequences from the CIF LMDB into a single FASTA (the input to
the sequence id map, clustering, and template search).

- Recipe: `structcooker.workflows.exports.fasta` (`fasta`)
- Config (extract): `db/metadata/cif_fasta.yaml` — `cif_pdb.lmdb` →
  `BIOMOL_ROOT/materials/intermediate/fasta/cif_pdb.fasta`
  (`reader.deserializer: ...codecs.from_bytes`, `writer.materializer: ...projections.write_fasta`)

```bash
pixi run structcooker build metadata/cif_fasta
```

The same recipe extracts the other sets' FASTAs: `metadata/teddymer_fasta`,
`metadata/ofd_long_fasta`, `metadata/ofd_short_fasta`,
`metadata/afm_homodimer_fasta`, `metadata/afm_heterodimer_fasta` (the AFM
halves are joined by `metadata/afm_fasta`, `metadata.concat_fasta`).
Train/valid variants use recipes matching their record layout:
`exports.fasta_attached` for `metadata/train_fasta` (attached train DB) and
`exports.fasta_transformed` for `metadata/valid_1_fasta` (the fanned-out valid-1
DB). Both feed valid stage 2. The older `exports.tv_fasta` has no `db/` config.

## Interacting sequences

- `exports.interacting_seq_ids` (`filtered_seq_ids`) /
  `db/metadata/interacting_seq_ids.yaml` (extract) — seq id pairs found in
  interfaces of `cif_pdb_attached.lmdb` → `metadata/interacting_seq_ids.tsv`.
- `exports.interacting_seq_clusters` (`interacting_seq_clusters`) /
  `db/metadata/interacting_seq_clusters.yaml` (materialize) — maps those ids to
  their `seq_cluster30_PDB.tsv` clusters → `metadata/interacting_seq_clusters.tsv`
  (a valid stage 2 input).

```bash
pixi run structcooker build metadata/interacting_seq_ids
pixi run structcooker build metadata/interacting_seq_clusters
```

## Per-chain CIF DB

`exports.cif_chain` (`biomoldict`) / `db/chain/cif_chain.yaml` (schema `C`,
rebuild) explodes each `cif_pdb.lmdb` record into one record per chain, keyed
`<pdbid>_<chain>` (best-occupancy assembly; the optional
`reference_chain_selection.lmdb` preserves legacy ties). All template builds read
template coordinates from it.

```bash
pixi run structcooker build chain/cif_chain
```

## Lists for template builds

Small materialize ops that write the work lists of the
[template ingest](ingest-template.md):
`exports.pdb_polypeptide_l` / `db/metadata/pdb_polypeptide_L.yaml`,
`exports.template_chain_list` / `db/metadata/template_chain_filelist.yaml`,
`exports.seqid_template_list` / `db/metadata/seqid_template_filelist.yaml`.

## Not in the release build

`exports.edge_node` (`monomer_map`, `interface_map`; monomer/interface graph
tables for graph splitting and profiling) and `exports.ccd` (`each_cif_lines`;
re-emits component CIF lines from the CCD LMDB) remain in the package but have
no `db/` config.
