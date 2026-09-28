# Verified PDB builds

The [5M storage policy](lmdb-storage-policy.md) governs finalization: larger DBs
publish immutable shard collections. Verified readers, indexing and content
fingerprints support both formats. Existing completed PDB outputs are retained.

Use `structcooker pdb-build` for the PDB CIF and supplied-MSA dependency graph in
`db/MANIFEST_pdb.yaml`. This entry point adds durable state and mandatory validation
to native DataCooker SLURM jobs. `build` and `build-all` retain their older semantics;
their existing-output checks do not provide these guarantees.

## Submit and resume

Install the locked Pixi environment and provide the inputs described in
[data provenance](data-provenance.md), including sequence-ID/cluster references,
chain selection and SignalP results. Supplied distillation sequences may be needed
as reference inputs to the shared sequence-ID space; this command does not build
the distillation datasets or generate missing MSA searches.

From the repository root:

```bash
export DATA_ROOT=/path/to/inputs
export OUTPUT_ROOT=/path/to/new/pdb-release
export PDB_RUN_DIR=/path/to/durable/pdb-run
export DATACOOKER_MAX_NODES=18
export DATACOOKER_EXCLUDE_NODES=node20
sbatch dev/pdb_build.sbatch

# Lightweight status; safe on the login node:
pixi run structcooker pdb-status --run-dir "$PDB_RUN_DIR"
```

The submission job plans stages and verifies outputs on its compute allocation;
native arrays, merge and index jobs use SLURM dependencies. Up to two independent
stages can progress concurrently. `DATACOOKER_MAX_NODES` limits each native pipeline,
not the combined allocation; SLURM QOS limits still apply. Worker defaults are
112 CPUs and 490 GB per node; override `DATACOOKER_CPUS_PER_NODE` and
`DATACOOKER_NODE_MEM_GB` when appropriate. The wrapper excludes node20 so it remains
available for experiments. Do not launch a second full build beside an existing
full-capacity run merely to validate the software.

Resubmit the same command with the same roots, run directory, inputs and software
to resume. Completed stages are reused only after their input and output content
identities match. Interrupted waits reuse recorded job IDs. An interruption during
validation resumes validation. Dependency stages are admitted only after their
parents pass verification. A nonzero exit is not a completed build.

## Completion contract

- Resolved configs, source packages and installed-package versions identify a run.
  Submitted jobs import a frozen source bundle. Input contents are SHA256 hashed;
  changed inputs, config, software or exception policy invalidate reuse.
- Inputs are hashed again after verification. Mutation during a build prevents a
  completion receipt. Hashing is deliberately full-content I/O; there is no unsafe
  size/mtime shortcut. Keep input snapshots immutable for the duration of a run.
- Every LMDB record is decoded and schema checked. Index keys, record sizes and
  metadata totals must agree. CIF checks cover every assembly's coordinate shape,
  chain references and bond endpoints. MSA checks cover array/header dimensions,
  value ranges, profiles and deletion means recomputed from the stored alignment.
- Native shard ledgers account for attempts, written records, empty filters and
  failures. Raw ingest keys must match the planned inputs after approved exclusions.
  A filtered rebuild may intentionally remove records; its write/skip counts are
  accounted for. Scientific selection correctness also relies on recipe tests.
- Only the five exact raw CIF inputs in `db/pdb-exceptions.json` may fail. Each
  waiver includes a source-content SHA256 and explanation. New failures stop the
  build. Do not add exceptions solely to obtain a green run.
- TSV/FASTA projections currently require declared, nonempty outputs and content
  hashes. They do not yet receive a complete independent semantic audit.

## Failure and ownership

`RUN_DIR/nodes/<config>/state.json` contains the phase, terminal job ID, attempt,
error and verification result. Attempts retain resolved engine configs, submission
logs, native scripts and shard reports. `summary.json` reports graph progress.

Normal and failed native script exits persist `JOB_ID.exit.json`; this permits
recovery after SLURM expires its job record. A killed job without an exit receipt
and without a scheduler record remains unconfirmed and cannot be silently reused.

An ambiguous `submitting` state requires reconciling the saved submission logs and
scheduler jobs before changing state. The command deliberately does not guess or
submit duplicates. A confirmed terminal failure permits a new attempt; unknown
status does not. Resolve errors before resubmitting.

Per-output locks prevent concurrent drivers. Persistent `.build-owner.json` files
also prevent a new run directory from bypassing ownership after the original
driver exits while native jobs remain alive. Resume the owner run. Do not remove
ownership/state files to force a retry with outstanding jobs.

Existing outputs without receipts are rejected. Prefer a fresh `OUTPUT_ROOT`.
`--rebuild-untracked` explicitly permits rebuilding untracked outputs; it does not
import old outputs as verified and cannot override another run's ownership. Keep
the existing recovered dataset untouched. These are build outputs, not an atomic
live-release switch: consumers should use a separately selected completed root.

## Acceptance evidence and limits

The automated tests cover interrupted waits and verification, uncertain submission,
confirmed failure, content changes/corruption, source snapshots, exact failure
coverage, exit receipts and persistent output ownership. Run `sbatch dev/test.sbatch`
on this cluster, not on the login node.

`dev/pdb_acceptance.py --run-dir /new/acceptance/path` runs inside a compute allocation
and submits real isolated native jobs. It builds two actual PDB CIF entries, three
CCD fixtures and a small protein MSA with a depth-capped derivative. It interrupts
the CLI, confirms job-ID reuse, confirms unchanged-stage reuse, then changes only
the MSA input and checks downstream invalidation. It temporarily uses
`db/_production_acceptance`; run only one acceptance instance at a time.
Run the repository-wide tests after acceptance finishes, since their configuration
portability checks enumerate every YAML under `db/`.
The provided submission wrapper is `sbatch dev/pdb_acceptance.sbatch`.

This acceptance test is not a fresh full-scale rerun of the complete PDB snapshot.
The prior bulk-output evidence remains in [CIF recovery](recovery-2026-09-10.md) and
[MSA reproduction](msa-reproduction.md). In particular, missing MSA coverage and
RNA source-version differences described there remain explicit data limitations.
Template search, distillation completion, fresh-environment deployment and atomic
release publication are outside this PDB build acceptance scope.

The September 11 acceptance run is recorded in
`logs/production/acceptance_v2/report.json` (SLURM job 266883). Restart preserved
terminal jobs 266886 and 266889; unchanged nodes were reused; changing the MSA input
rebuilt only the MSA and its depth-capped derivative. Job 266900 passed all 167
tests after this run. Static checks passed for both source packages with no Pyright
errors or warnings. The final accounting check and regression run are recorded
separately under `logs/production/final_*` and
`logs/production/acceptance_v2/final-verification.json`.
