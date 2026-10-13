# Tutorials

These tutorials cover the **domain workflows**. For how the recipe engine and
LMDB builds work, see the
[DataCooker docs](https://CSSB-SNU.github.io/DataCooker/). The
[data flow](data-flow.md) page shows how they fit together end to end.

Every output is described by exactly one `db/<group>/<name>.yaml` (recipe,
reader/writer hooks, source, target, and a `schema:` tag) and built with

```bash
pixi run structcooker list                       # every db config: op, schema, target
pixi run structcooker build <group>/<name> --show     # resolved invocation, no submit
pixi run structcooker build <group>/<name> --dry-run  # plan + write sbatch, no submit
pixi run structcooker build <group>/<name>            # plan, size, submit (afterok-chained)
```

The op (`build` / `rebuild` LMDB pipeline, or the `materialize` / `extract` /
`parallel` workflow ops) is inferred from the config's shape, and the SLURM jobs
are sized and submitted by `structcooker` itself — there is no map size, shard
count, or sbatch script to write. Paths resolve from one root, `BIOMOL_ROOT` (default
`/data/shared/cssb_data/BioMol_clean`); only the shared search databases come from
`SEQ_DB_ROOT` (default `/data/shared/cssb_data/seq_db`).
The whole release is the DAG in `db/MANIFEST_all.yaml`; to build it end to end
see the [end-to-end build](e2e-build.md) runbook and, for where inputs live,
[raw materials](raw-materials.md).

## Ingest (build LMDBs)

- [CCD ingest](tutorials/ingest-ccd.md) — Chemical Component Dictionary.
- [CIF ingest](tutorials/ingest-cif.md) — RCSB mmCIF structures → `CIFMol`.
- [MSA ingest](tutorials/ingest-msa.md) — a3m alignments and depth caps.
- [Template ingest](tutorials/ingest-template.md) — template hits → `TemplateMol`.
- [OpenFold3 distillation ingest](tutorials/openfold-distillation.md) —
  MSA / structure / template from the distillation datasets.

## Build the inputs

- [Search](tutorials/search.md) — sequence clustering, MSA search, template search.
- [Metadata](tutorials/metadata.md) — sequence id map, attach metadata.

## Curate & inspect

- [Filters](tutorials/filters.md) — train / valid splits, MSA cleanup.
- [Exports](tutorials/exports.md) — FASTA, interacting sequences, per-chain DB.
- [Analysis](tutorials/analysis.md) — profiling and validation.
