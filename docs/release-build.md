# Verified native SLURM releases

`structcooker release-build` submits finite planning, worker, publication, index,
and verification jobs. Only a verified stage releases its dependents. The command
returns after submission; `summary.json` with `complete: true` records completion.
There is no resident controller. `build-all --workdir RUN` uses this same path.

For a new distillation release, supply existing CCD/chain/sequence reference
artifacts separately from a fresh destination:

```bash
export DATA_ROOT=/path/to/raw-data
export REFERENCE_ROOT=/path/to/reference-release
export OUTPUT_ROOT=/path/to/new-release
export DATACOOKER_EXCLUDE_NODES=node20
pixi run --frozen structcooker release-build --run-dir /path/to/run-state
```

Inspect input availability first with `structcooker preflight --help` and the
manifest-specific options. `db/MANIFEST_distillation.yaml` declares all configured
CIF, attached CIF, contacts, MSA depth variants and template stages. It assumes the
reference artifacts already exist; it does not download the raw archive or create
CCD/chain reference databases. For PDB use `--manifest db/MANIFEST_pdb.yaml
--policy db/pdb-exceptions.json` with the same command.

Repeat the exact command with the same roots, code, lock and run directory to
resume. Active recorded jobs are retained. Terminal attempts reuse worker shards
and retry publication/indexing; a completed run rechecks content identities on
compute nodes before reuse. Changed definitions require a new run and output
root. An ambiguous submission or unknown scheduler state fails closed and needs
inspection of the recorded job manifests; it is not automatically resubmitted.
Do not delete run-state or physical shards behind a published shard database.

Reader maps and recovery overrides are explicit configuration inputs. Each
attempt pins its reader references; input and output hashes protect reuse.
Duplicate raw keys fail by default. Distillation MSA configs explicitly choose the
lexicographically first source path and retain an `input-selection.json` ledger.
A failed selected source remains a failure; another duplicate is not silently
substituted. Approved raw failures require an exact key, SHA256 and reason in the
policy JSON. Template records retain missing-chain and selection ledgers;
decode/transform errors fail the record.

Outputs with at most 5,000,000 records are merged; larger outputs use `shard_lmdb`.
`pixi.lock` is part of the deliverable and scheduler jobs use `pixi run --frozen`.
Keep it with both source packages. The older `pdb-build` synchronous interface
remains for compatibility; the release interface above is the finite SLURM path.

Validation on 2026-09-17: a two-input MSA fixture completed base → depth cap →
verification on SLURM, then reused the same attempts on repeat invocation. Removing
the fixture cap index and resuming restored a verified output from the same shards.
This checks the scheduler and recovery path, not full-scale capacity or fresh-machine
installation. Existing production jobs continue using their previously frozen code.
