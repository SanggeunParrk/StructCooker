# StructCooker

StructCooker is a domain package built on top of DataCooker for structural biology workflows.

## Package Layout

- `src/structcooker/models`: domain objects and structural containers
- `src/structcooker/readers`: file and LMDB readers
- `src/structcooker/ops`: reusable transformation logic
- `src/structcooker/workflows`: DataCooker recipe modules
- `src/structcooker/writers`: output materializers

## Repository Layout

- `configs/ingest`: LMDB build configs
- `configs/filters`: filtering and cleanup workflows
- `configs/metadata`: metadata enrichment workflows
- `configs/exports`: export and extraction workflows
- `configs/search`: MSA and template search workflows
- `configs/analysis`: database inspection and validation workflows
- `submits/*`: cluster entrypoints grouped by the same domains
- `submits/legacy`: environment-specific legacy scripts kept for reference
- `scripts/dev`: ad-hoc inspection helpers
- `scripts/maintenance`: one-off maintenance utilities
