# Analysis

Analysis workflows validate and profile the databases produced by the pipeline.

## Release checks

Two CLI checks cover a built release; both submit SLURM jobs.

```bash
# Is the release tree complete and self-contained? Every stage's output at its
# configured path, shards inside their own collection, entry counts matching
# .meta.json / .index.tsv, and a sample of records decoding cleanly.
pixi run structcooker release-check            # --manifest db/MANIFEST_all.yaml (default)
# -> OUTPUT_ROOT/release_check.json

# Does a from-scratch build reproduce the release? Build the whole DAG on a
# small sample of inputs and compare every sampled output with the release.
pixi run structcooker e2e-sample prepare <ROOT>
pixi run structcooker e2e-sample run <ROOT>
pixi run structcooker e2e-sample status <ROOT>
pixi run structcooker e2e-sample verify <ROOT>  # -> <ROOT>/verify.json
```

See the [end-to-end build](../e2e-build.md#5-sampled-end-to-end-check) for the
sample sizes and how to resume. Before building, `structcooker inspect --strict
--manifest db/MANIFEST_all.yaml` checks that every external input is present
and the DAG is closed.

## Recipe-level checks

The recipes below remain in `structcooker.workflows.analysis` but have no `db/`
config and are not part of the release build; run them from Python.

- `analysis.db_profile` (`monomer_clusters`, `interface_clusters`) summarises
  the cluster composition of a built database.
- `analysis.template_check` (`results`) reads a template LMDB and reports
  per-entry sanity (resolved chains, alignment coverage).
- `analysis.ccd_validate` (`issues`) runs per-component checks over the CCD
  LMDB:

```python
from datacooker import execute
from structcooker.workflows.analysis.ccd_validate import RECIPE, TARGETS

issues = execute(RECIPE, {"chem_comp": entry}, targets=TARGETS)["issues"]
```
