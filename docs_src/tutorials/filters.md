# Filters

Filtering workflows derive curated subsets (train / valid splits, cleaned MSAs)
from the core databases. Structure filters use the **rebuild** path
(`old_env_path` → `new_env_path`) with cutoffs in the config's `parameters`.

## Train / valid structure filtering

`structcooker.workflows.filters.data` (`cifmol_dict`) keeps structures passing
resolution / date / size cutoffs.

| Config | Window | Cutoffs (example) |
| --- | --- | --- |
| `configs/filters/train_filter.yaml` | `end_date: 2021-09-30` | resolution ≤ 9.0, ≤ 10240 tokens, ≤ 300 chains |
| `configs/filters/valid_filter_1.yaml` | 2021-09-30 → 2023-09-30 | resolution ≤ 4.5, ≤ 2560 tokens, ≤ 20 chains |

```bash
sbatch submits/filters/train_filter.sh
sbatch submits/filters/valid_filter_1.sh
# directly:
pixi run python -m datacooker.cli.lmdb rebuild configs/filters/valid_filter_1.yaml
```

Stage-2 validation (`filters.validation_stage1` / `validation_stage2` /
`validation_stage2_metadata`) further refines the valid split.

## MSA filtering

- `filters.a3m` (`sequences`, `headers`) / `configs/filters/a3m_filter.yaml`
  — depth/quality filter on alignments.
- `filters.a3m_cleanup` (`results`) /
  `configs/filters/remove_lower_from_a3m.yaml` — drop insertion (lowercase)
  columns.

```bash
pixi run python -m datacooker.cli.lmdb rebuild configs/filters/a3m_filter.yaml
```

## Disordered subset

`filters.disordered` (`cifmol_dict`) / `configs/filters/disordered_data_filter.yaml`
selects the disordered-region subset.
