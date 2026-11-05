# abag_distillation

The antibody-antigen self-distillation set built from **ASD v2.0, AbRank and POISK**: pair table `abag/pairs/pairs_basic.tsv`, antigens with an MSA and templates in two versions (released up to 2021-09-30 and up to 2026),
and the final table `abag/final/pairs_final.tsv` with a seq_id and a cluster for every chain.

| File | What |
|---|---|
| [overview.html](overview.html) | one page: what was done and how (flow, funnel, evidence, decisions); open in a browser |
| [RECIPE.md](RECIPE.md) | the recipe: stages in order, commands, inputs and outputs, how to change it |
| [run.sh](run.sh) | runs the stages (`ids`, `pairs`, `antigen`, `final`, `all`) |
| [sources.md](sources.md) | the three sources and their helpers, sequences, ids, antigen MSA and templates, CDR clustering |
| [affinity-filter.md](affinity-filter.md) | the Kd <= 1 uM rule per source and `affinity_type`, what is left out as ambiguous |
| [light-chain-pairing.md](light-chain-pairing.md) | SelfDistillation's light-chain recovery, our SAbDab measurement, what we do |
| [pair-filters.md](pair-filters.md) | duplicates, chain status, polyprotein, experimental structure, the antigen filter, with counts |

The stage configs are in `db_configs/abag/`, the recipes in `src/structcooker/workflows/metadata/abag_*.py`, the readers in `src/structcooker/instructions/readers/abag*.py` (and `asd.py`, `asd_binding.py` for ASD).
