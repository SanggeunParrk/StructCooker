# StructCooker

StructCooker is a domain package built on top of DataCooker for structural biology workflows.

## Package Layout

- `src/structcooker/mols`: BioMol-backed structural domain objects
- `src/structcooker/utils`: shared mapping and domain utilities
- `src/structcooker/instructions/readers`: external data loaders and LMDB/file readers
- `src/structcooker/instructions/transforms`: adapters, serializers, and reusable transformation instructions
- `src/structcooker/instructions/writers`: output materializers

## Workflow Layout

- `workflows/ingest`: LMDB build recipes
- `workflows/filters`: dataset filtering and cleanup recipes
- `workflows/metadata`: metadata enrichment and export-prep recipes
- `workflows/exports`: downstream artifact extraction recipes
- `workflows/search`: MSA, template, and clustering recipes
- `workflows/analysis`: validation and profiling recipes
