# Analysis

Analysis workflows validate and profile the databases produced by the pipeline.

## Database profiling

`structcooker.workflows.analysis.db_profile`
(`monomer_clusters`, `interface_clusters`) summarises cluster composition of a
built database.

- Configs: `configs/analysis/analyze_train.yaml`, `analyze_valid1.yaml`, …

```bash
sbatch submits/analysis/analyze_db.sh
```

## Template database check

`structcooker.workflows.analysis.template_check` (`results`) reads the template
LMDB and reports per-entry sanity (resolved chains, alignment coverage).

- Config: `configs/analysis/check_template_db.yaml`

```bash
pixi run python -m datacooker.cli.workflow extract-lmdb configs/analysis/check_template_db.yaml
```

## CCD validation

`structcooker.workflows.analysis.ccd_validate` (`issues`) runs per-component
checks over the CCD LMDB:

```python
from datacooker import execute
from structcooker.workflows.analysis.ccd_validate import RECIPE, TARGETS

issues = execute(RECIPE, {"chem_comp": entry}, targets=TARGETS)["issues"]
```
