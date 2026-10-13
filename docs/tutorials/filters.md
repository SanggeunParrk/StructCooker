# Filters

Filtering workflows derive curated subsets (train / valid splits, cleaned MSAs)
from the core databases. Structure filters are `rebuild` ops
(`old_env_path` → `new_env_path`) with cutoffs in the config's `parameters`.

## Train / valid structure filtering

The filter recipes trim each assembly by token count, polymer chain count,
resolution and date window, remove water, and drop SignalP-flagged sequence
regions (`metadata_recipe: metadata.load_signalp`, reading
`materials/reference/signalp` and `metadata/seq_id_map.tsv`).

- `structcooker.workflows.filters.data_attached` (`cifmol_attached_dict`) — train
  splits, read from the **attached** DB so records keep their seq and cluster ids.
- `structcooker.workflows.filters.data` (`cifmol_dict`) — valid stage 1, read
  from the base `cif_pdb.lmdb`.

| Config | Source | Window | Cutoffs |
| --- | --- | --- | --- |
| `db/train/train_20210930.yaml` | `cif_pdb_attached.lmdb` | `end_date: 2021-09-30` | resolution ≤ 9.0, ≥ 1 token (no max), ≤ 300 polymer chains |
| `db/train/train_20260301.yaml` | `cif_pdb_attached.lmdb` | `end_date: 2026-03-01` | resolution ≤ 9.0, ≥ 1 token (no max), ≤ 300 polymer chains |
| `db/valid/valid1.yaml` (schema `I`) | `cif_pdb.lmdb` | 2021-09-30 → 2023-09-30 | resolution ≤ 4.5, ≤ 2560 tokens, ≤ 20 polymer chains |

```bash
pixi run structcooker build train/train_20210930
pixi run structcooker build valid/valid1
```

Validation then continues in two steps: `db/valid/valid1_attach.yaml` attaches
sequence metadata to the valid-1 DB (see [Metadata](metadata.md)), and
`db/valid/valid2.yaml` (`filters.validation_stage2` /
`validation_stage2_metadata`, `cifmol_attached_dict`) drops entries whose
clusters overlap the train or valid-1 sequences (`metadata/train_fasta`,
`metadata/valid_1_fasta`) using `seq_cluster30_PDB.tsv` and the interacting
clusters from [Exports](exports.md#interacting-sequences).

```bash
pixi run structcooker build valid/valid1_attach
pixi run structcooker build valid/valid2
```

`db/train/train_20210930_mpnn.yaml` (`filters.mpnn`, schema `H`) is a
model-input view of the 2021-09-30 train DB: it keeps only the features the MPNN
and cropper read and explodes each `{assembly, model, altloc}` into its own key
(see [MPNN pruning](../mpnn_pruning.md)).

## MSA filtering

`filters.a3m_cleanup` (`results`) drops insertion (lowercase) columns from a3m
files — the hmmbuild input for template search. It runs as a `parallel` op:
`db/template/msa_wo_lower.yaml` (PDB), `db/teddymer/msa_wo_lower.yaml`,
`db/afdb_multimer/msa_wo_lower.yaml`.

```bash
pixi run structcooker build template/msa_wo_lower
```

MSA depth is capped by `exports.cap_msa_depth` rebuilds (see
[MSA ingest](ingest-msa.md#depth-caps)). The older depth/quality filter
`filters.a3m` has no `db/` config and is not part of the release build.

## Disordered subset

The disordered subset is now the OpenFold disordered distillation set
(`db/distillation/disordered_*.yaml`, see
[OpenFold3 distillation ingest](openfold-distillation.md)). The older
`filters.disordered` (`cifmol_dict`) recipe has no `db/` config and is not part
of the release build.
