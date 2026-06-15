# Metadata

Metadata workflows assign stable sequence ids and attach per-sequence metadata
onto structure records.

## Sequence id map

Assign a stable id to every unique sequence in the extracted FASTA (an existing
map can be reused so ids stay stable across rebuilds).

- Recipe: `structcooker.workflows.metadata.seq_id_map` (`seq_id_map`)
- Config: `configs/metadata/build_seq_id_map.yaml`
- Writer: `...writers.projections.write_seq_id_map` → `seq_id_map.tsv`

```bash
pixi run python -m datacooker.cli.workflow run configs/metadata/build_seq_id_map.yaml
```

## Attach sequence metadata

Re-build a CIF LMDB into an *attached* LMDB, enriching each `CIFMol` with
sequence-level metadata (ids, clusters, ...). This uses the **rebuild** path:
an `old_env_path` is read through a `reader.adapter`, a `metadata_recipe` loads
shared metadata once, and records are written to `new_env_path`.

- Recipe: `structcooker.workflows.metadata.attach` (`cifmol_attached_dict`)
- Metadata recipe: `metadata.load_sequence_metadata`
  (`seq_metadata_map`, `seqid2seq`, `seqclusters2seqids`)
- Config: `configs/metadata/attach_seq_metadata_train.yaml`

```bash
sbatch submits/metadata/attach_seq_metadata.sh
# directly:
pixi run python -m datacooker.cli.lmdb rebuild configs/metadata/attach_seq_metadata_train.yaml
```

## Extract metadata

`structcooker.workflows.metadata.extract` (`metadata_dict`) /
`configs/metadata/extract_metadata.yaml` dumps a metadata table for downstream
filtering and analysis.
