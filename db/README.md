# `db/` — one BioMol DB, one declarative config

Each `db/**/<name>.yaml` fully describes how to build **one** BioMol LMDB. It is a
[datacooker](../libs/datacooker) *engine config* (recipe / reader / writer + source +
target) plus a single **`schema:`** tag. Nothing else — no `n_jobs`, `chunk_size`,
`file_list`, or size tiers. Those are derived, not authored.

```
structcooker build <name>        # plan + build (op auto-inferred), afterok-chained
structcooker build <name> --dry-run   # plan + write sbatch, do not submit
structcooker build <name> --show      # print resolved schema/E + invocation, exit
structcooker list                # every db/*.yaml with op, schema, target
```

## How a build runs

`structcooker build` reads the `schema:` tag, looks up that schema's measured
**expansion factor E** (peak-live-mem ÷ stored-value-bytes) in
[`schemas.py`](../src/structcooker/schemas.py), then hands the config to
`datacooker pipeline`, which:

1. **sizes** every item — source-file `st_size` (build) or stored value bytes from
   `<db>.index.tsv` (rebuild);
2. **plans** memory-safe size tiers from those sizes and E (`peak = n_jobs × E ×
   max_bytes`), so giants never cluster onto a worker and OOM;
3. **submits** one job array per tier → **merges** the shards → writes the target's
   own `.index.tsv` + `.meta.json`, all `afterok`-chained.

No human tunes n_jobs / mem / chunk / shards. The `schema:` tag never reaches the
engine: `structcooker` materializes a schema-stripped `engine.yaml` in the workdir.

## Op is inferred

| config has            | op        | source sizing            |
|-----------------------|-----------|--------------------------|
| `env_path` + `data_dir`/`file_pattern` | `build`   | source-file `st_size`    |
| `old_env_path` + `new_env_path`        | `rebuild` | `<source>.index.tsv`     |

## Schemas (see `schemas.py`)

| schema | value                              | key           | E     |
|--------|------------------------------------|---------------|-------|
| A      | `{assembly_dict, metadata_dict}`   | pdbid         | 100   |
| B      | cifmol_attached (per assembly)     | pdbid         | 100   |
| C      | chain BioMol                       | pdbid_chain   | 100   |
| D      | `{template_mols}`                  | pdbid_chain   | 20    |
| E      | `{msa_dict}`                       | seqid         | 40    |
| F      | raw scalar bytes                   | mixed         | 1     |

## Targets

Every DB writes under **`/data/shared/cssb_data/BioMol_test/`** — the original
`BioMol/` tree is left untouched until a rebuilt DB is verified against it.

## Layout

```
db/
  pdb/     cif.yaml  cif_attached.yaml          # A ingest, B attach
  train/   train_20210930.yaml  train_20260301.yaml   # B release-date filters
  chain/   …   template/   …   msa/   …   ccd/   …
  distillation/   …        # code/config only — build deferred
```
