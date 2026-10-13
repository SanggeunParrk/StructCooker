# Development and verification

Run from the repository root after installing the Pixi environment and initializing
`libs/datacooker`. Static checks do not need a SLURM allocation:

```bash
pixi run ruff check src tests libs/datacooker/src libs/datacooker/tests
pyright -p pyproject.toml
```

Pyright is a separately installed development tool; the repository configuration
resolves dependencies from `.pixi/envs/default` and checks both packages' source.
It no longer depends on a temporary recovery configuration.

The full suite includes native CCD ingestion and Ray execution. Run it on a compute
node, using the partition and QoS for your cluster. For this deployment:

```bash
mkdir -p logs/tests
sbatch -p cpu-short --qos=cpu-short-q -c 8 --mem=32g -o logs/tests/structcooker_%j.out \
  --wrap "pixi run --frozen python -m pytest -q -p no:cacheprovider tests"
# DataCooker's own suite, from its directory (its Ray workers import `tests.*`):
sbatch -p cpu-short --qos=cpu-short-q -c 8 --mem=32g -o logs/tests/datacooker_%j.out \
  --wrap "cd libs/datacooker && pixi run --frozen --manifest-path ../../pyproject.toml python -m pytest -q -p no:cacheprovider tests"
```

Both use temporary test data and never touch a release tree. Check each job's terminal
state and the pytest summary line; a successful submission alone is not a passed test.

To check a dataset configuration without submitting any builds:

```bash
pixi run structcooker inspect --manifest db/MANIFEST_cifcore.yaml --strict
```

This checks nested/list inputs, CCD inputs, supplied model choices, duplicate output
paths, cycles, unknown dependencies, and missing producer dependencies. An existing
path is not proof of valid contents, and the tools inventory is informational.

## Repository artifacts

- `src/`, `db/`, and `libs/datacooker/src/` are the supported implementation and configs.
- `tests/fixtures/` contains small distributable test inputs; tests need no recovery snapshot.
- `logs/recovery/` preserves the September 10 audit evidence, frozen source snapshots,
  and one-time recovery scripts. Normal builds and tests do not import from it.
- `scripts/` holds only `scripts/staging/` (raw-input staging, see the
  [end-to-end build](e2e-build.md)). The legacy `configs/` and `submits/`, the
  one-off `scripts/{dev,maintenance}`, and the `dev/` and `scratchpad/` scratch
  areas were archived out of the checkout on 2026-10-13; use `db/` and the
  `structcooker` CLI for builds.
- `BioMol_clean/repair_backups/` contains recovery backups. The supplied CCD and
  model-choice references under `BioMol_clean/metadata/` are reproduction inputs,
  not disposable caches.

Changes remain uncommitted in both the root repository and the DataCooker submodule.
