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
test. Tests also keep `db_configs/MANIFEST_all.yaml` equal to the union of the per-set manifests and
every `db_configs/` config resolvable.

To check configs without submitting anything:

```bash
pixi run structcooker inspect --strict --manifest db_configs/MANIFEST_all.yaml
pixi run structcooker build-all --manifest db_configs/MANIFEST_all.yaml --dry-run
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
db_configs/                         one config per DB, in a folder per set: db_configs/<group>/<name>.yaml
db_configs/MANIFEST_*.yaml          the release DAG (per set) and MANIFEST_all.yaml (the whole release)
db_configs/release-exceptions.json  inputs a release stage may fail on: key, reason, SHA256
db_configs/<group>/params/          parameter tables of a DB (not engine configs; `structcooker list` skips them), e.g. interface/params/
src/structcooker/
  cli.py                    the `structcooker` commands
  paths.py                  BIOMOL_ROOT / SEQ_DB_ROOT
  build/                    build machinery: build_state (fingerprints, output claims, locks), build_verify (full output
                            verification), preflight (input/output paths, dependency errors), workflow_ops (op inference and the
                            SLURM job), schemas (value type and expansion factor), downloads, cluster (node size)
  release/                  release level: dag (build-all: frozen definition, plan / verify jobs per stage, resume), check
                            (release-check), docs (release-docs, pages from docs/release/), manifests (MANIFEST_all.yaml = the
                            union of the per-set manifests), key_contracts, sample_e2e (e2e-sample), foldbench
  instructions/readers/     loaders and key builders: cif, a3m, openfold, teddymer, afdb, lmdb, asd, abag
  instructions/transforms/  adapters, codecs, msa, template, pruning, sequence, ...
  instructions/writers/     projections: FASTA / TSV materializers
  instructions/mols/        CIFMol, CIFMolAttached, TemplateMol
  instructions/interface/   the interface-unit DB (detect, duplex, adapt, ...; its parameters are db_configs/interface/params/)
  instructions/design/, fragdb/, potential/   the design set, the fragment DB, the potential frames and pairs
  workflows/                recipes: ingest/, metadata/, search/, filters/, exports/
  utils/                    seq_id and DB-code helpers
libs/datacooker/            the engine (submodule)
scripts/cluster/cssb.env    this cluster's settings
scripts/staging/            raw-input staging (e2e-build §2), build_inference_ccd.py, consolidate_shards.py
docs/release/               hand-written pages of the released data tree and scenarios.yaml (release-docs copies them)
```

The engine has two layers, **instructions** (the functions) and **workflows** (the recipes a config points at with `recipe_path`); everything else in `src/structcooker/` is the machinery
around them. A config is `db_configs/<group>/<name>.yaml`; nothing about a DB is configured anywhere else.

`logs/` keeps job logs and recovery evidence; nothing in a build imports from it.

## Quality gates

Every push runs the same four checks (`.github/workflows/ci.yml`); run them before you commit:

```bash
pixi run ruff check src tests scripts          # lint: select = ALL minus the principled ignores in pyproject.toml
pixi run ruff format --check src tests scripts
pixi run pyright                               # types of src/ (the engine is checked in its own repository)
pixi run pytest tests -q --cov=src/structcooker   # 812 tests, 66% of the statements (the recipes in workflows/ are glue, covered by e2e-sample)
```

- **Library code (`src/`)** has no suppressions beyond a justified `# noqa: <rule>  (reason)` on the line (subprocess with a fixed argv, seeded sampling that is not security, a signature the engine fixes).
- **Tests** may use bare asserts, skip docstrings and annotations, touch private helpers (`SLF001`) and use fakes with the signature of what they replace (`ARG001/2`).
- **Scripts (`scripts/`)** are one-off drivers, not library code: no annotations or docstrings, `print`, `open()` and short names are fine; correctness and security rules stay on.
- The recipes (`workflows/`) are glue the engine runs from a config; `tests/workflows/test_recipes.py` checks that each loads, declares `RECIPE` / `TARGETS`, and is used by a config. Their behaviour is
  covered by `structcooker e2e-sample` and `release-check`, not by the unit-test coverage figure.

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
