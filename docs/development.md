# Development

## Checks and tests

Static checks need no allocation:

```bash
pixi run ruff check src tests libs/datacooker/src libs/datacooker/tests
pyright -p pyproject.toml          # installed separately; resolves from .pixi/envs/default
```

The test suites include native CCD ingestion and Ray execution, so run them as batch jobs,
never on a login node. With this cluster's short partition:

```bash
mkdir -p logs/tests
sbatch -p cpu-short --qos=cpu-short-q -c 8 --mem=32g -o logs/tests/structcooker_%j.out \
  --wrap "pixi run --frozen python -m pytest -q -p no:cacheprovider tests"
# DataCooker's own suite, from its directory (its Ray workers import `tests.*`):
sbatch -p cpu-short --qos=cpu-short-q -c 8 --mem=32g -o logs/tests/datacooker_%j.out \
  --wrap "cd libs/datacooker && pixi run --frozen --manifest-path ../../pyproject.toml python -m pytest -q -p no:cacheprovider tests"
```

Both use temporary data under `tests/fixtures/` and never touch a release tree. Check each
job's terminal state and the pytest summary line; a successful submission is not a passed
test. Tests also keep `db/MANIFEST_all.yaml` equal to the union of the per-set manifests and
every `db/` config resolvable.

To check configs without submitting anything:

```bash
pixi run structcooker inspect --strict --manifest db/MANIFEST_all.yaml
pixi run structcooker build-all --manifest db/MANIFEST_all.yaml --dry-run
```

`inspect` checks nested and list inputs, CCD inputs, supplied model choices, duplicate output
paths, cycles, unknown dependencies and missing producers. An existing path is not proof of
valid contents.

## Cluster settings

Nothing in the code names a partition, a node or a site path except as a default. The CSSB
values are in `scripts/cluster/cssb.env`; `source` it before building (jobs inherit the
submitting shell's environment), or copy it for another cluster.

| Variable | CSSB value | Used for |
|---|---|---|
| `BIOMOL_ROOT` | `/data/shared/cssb_data/BioMol_clean` | the release tree: `materials/`, `lmdb/`, `metadata/`, `docs/` |
| `SEQ_DB_ROOT` | `/data/shared/cssb_data/seq_db` | UniRef30 / MMseqs2 (only `teddymer/msa_search`, the AFM missing-MSA job) |
| `DATACOOKER_PARTITION`, `DATACOOKER_QOS` | `cpu-long`, `cpu-long-q` | every submitted job |
| `DATACOOKER_EXCLUDE_NODES` | `node01,node04,node07` | nodes whose `/dev/shm` cannot hold joblib's semaphores |
| `DATACOOKER_CPUS_PER_NODE`, `DATACOOKER_NODE_MEM_GB` | `112`, `490` | a full-node job (planner, parallel searches, verification) |

`STRUCTCOOKER_MAX_NODES` caps every config's `node_count` (the sampled e2e run sets 1).

## Code layout

```
db/                         one config per stage; MANIFEST_*.yaml (per set) and MANIFEST_all.yaml (release)
db/release-exceptions.json  inputs a release stage may fail on: key, reason, SHA256
src/structcooker/
  cli.py                    the `structcooker` commands
  release.py                build-all: frozen definition, plan / verify jobs per stage, resume
  build_verify.py           full output verification before a stage releases its dependents
  build_state.py            fingerprints, output claims, locks
  preflight.py              inspect: input and output paths of a config, dependency errors
  workflow_ops.py           op inference and the SLURM job for workflow ops
  release_check.py          release-check
  release_docs.py           release-docs (pages from release_docs/)
  sample_e2e.py             e2e-sample
  foldbench.py              foldbench
  schemas.py                schema tags: value type and expansion factor
  paths.py, cluster.py      BIOMOL_ROOT / SEQ_DB_ROOT, node size from the environment
  manifests.py              MANIFEST_all.yaml = union of the per-set manifests
  downloads.py              structcooker download
  mols/                     CIFMol, CIFMolAttached, TemplateMol
  instructions/readers/     loaders and key builders: cif, a3m, openfold, teddymer, afdb, lmdb
  instructions/transforms/  adapters, codecs, msa, template, pruning, ...
  instructions/writers/     projections: FASTA / TSV materializers
  workflows/                recipes: ingest/, metadata/, search/, filters/, exports/
libs/datacooker/            the engine (submodule)
scripts/cluster/cssb.env    this cluster's settings
scripts/staging/            raw-input staging (e2e-build §2), build_inference_ccd.py, consolidate_shards.py
release_docs/               hand-written release pages and scenarios.yaml
```

`logs/` keeps job logs and recovery evidence; nothing in a build imports from it.

## LMDB storage

- **Merge or shard.** After a build, the output records are counted. At most 5,000,000
  (`datacooker.lmdb.sharded.MERGE_LIMIT`) are merged into one LMDB; more are published as a
  shard collection: `<name>.lmdb/` holds `shards.json` (format `shard_lmdb`, relative shard
  paths, counts) instead of `data.mdb`. Publication checks that keys are unique across shards
  and never copies payloads. The count is of output records, so an exploding rebuild counts
  its exploded keys.
- **Read either with one call:** `datacooker.lmdb.sharded.open_env(path, readonly=True,
  lock=False)` (get, entry count, sorted forward iteration). Code that calls `lmdb.open`
  directly does not see a collection.
- **Sidecars.** Every DB has `<name>.lmdb.index.tsv` (key and stored size per record, what a
  rebuild is planned from) and `<name>.lmdb.meta.json` (counts). `release-check` requires
  them to agree with the DB.
- **Shards are permanent.** Keep shards immutable and inside their collection; move a
  collection as a whole. A patch adds a shard rather than rewriting one.
- **MSA depth variants** are separate DBs rebuilt from the base MSA DB (`*_d16k`, `_d5k`,
  `_d2k`, `_d512`), keyed identically: verification requires every source key in the cap.
- `map_size` defaults to 300 GB and is a ceiling, not an allocation; a config sets more where
  a merged DB is larger (`afdb_multimer/msa`: 1.5 TB).
