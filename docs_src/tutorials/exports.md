# Exports

Export workflows pull artifacts out of the LMDBs. Most read an LMDB through a
recipe (`cli.workflow extract-lmdb`); a few derive tables from existing files
(`cli.workflow run`).

## FASTA

Extract all sequences from the CIF LMDB into a single FASTA (the input to the
sequence id map, clustering, and MSA search).

- Recipe: `structcooker.workflows.exports.fasta` (`fasta`)
- Config: `configs/exports/extract_fasta_whole.yaml`
  (`reader.deserializer: ...codecs.from_bytes`, `writer.materializer: ...write_fasta`)

```bash
sbatch submits/exports/extract_fasta.sh
# directly:
pixi run python -m datacooker.cli.workflow extract-lmdb configs/exports/extract_fasta_whole.yaml
```

Train/valid-specific variants: `exports.tv_fasta` (`fasta`) /
`configs/exports/extract_fasta_train.yaml`, `extract_fasta_valid_1.yaml`.

## Interacting sequences

- `exports.interacting_seq_ids` (`filtered_seq_ids`) /
  `configs/exports/extract_interacting_seq_ids.yaml` — sequence ids found in
  interfaces.
- `exports.interacting_seq_clusters` (`interacting_seq_clusters`) /
  `configs/exports/extract_interacting_seq_clusters.yaml` — maps those ids to
  their clusters.

```bash
pixi run python -m datacooker.cli.workflow run configs/exports/extract_interacting_seq_clusters.yaml
```

## Edge / node tables

`exports.edge_node` (`monomer_map`, `interface_map`) /
`configs/exports/extract_edge_node_*.yaml` produce the monomer/interface graph
tables used for train/valid graph splitting and profiling.

## CCD CIF round-trip

`exports.ccd` (`each_cif_lines`) / `configs/exports/extract_ccd_cif.yaml`
re-emits component CIF lines from the CCD LMDB.
