# StructCooker

StructCooker is a domain package built on top of DataCooker for structural biology workflows.

## Layout

- `src/structcooker/models`: domain objects
- `src/structcooker/readers`: input adapters and LMDB readers
- `src/structcooker/ops`: workflow operations
- `src/structcooker/workflows`: DataCooker recipe modules
- `src/structcooker/writers`: output materializers
- `configs/`: runnable workflow configurations
- `submits/`: cluster/job entrypoints
