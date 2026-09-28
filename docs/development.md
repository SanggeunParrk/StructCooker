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
mkdir -p logs/quality
sbatch --partition=cpu-long --qos=cpu-long-q --exclude=node04 \
  --output=logs/quality/tests_%j.out --error=logs/quality/tests_%j.err dev/test.sbatch
```

`dev/test.sbatch` uses one node, starts an isolated local Ray instance, imports the
workspace packages, and propagates pytest's exit status. It uses temporary test data
and does not modify the production or recovered databases. Confirm the job's terminal
state and exit code; a successful submission alone is not a passed test.

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
- `dev/` contains maintained development entry points.
- `logs/recovery/` preserves the September 10 audit evidence, frozen source snapshots,
  and one-time recovery scripts. Normal builds and tests do not import from it.
- Legacy `scripts/`, `configs/`, and the inherited `scratchpad/` are retained for
  provenance and local work; use `db/` and the CLI for current builds.
- `BioMol_clean/repair_backups/` contains recovery backups. The supplied CCD and
  model-choice references under `BioMol_clean/metadata/` are reproduction inputs,
  not disposable caches.

Changes remain uncommitted in both the root repository and the DataCooker submodule.
