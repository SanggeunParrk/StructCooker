# Full distillation CIF build

**September 13 update:** the long merge DAG described below was replaced by the
[5M shard storage policy](lmdb-storage-policy.md). Old jobs 266824–266825 and
266848–266853 were cancelled; jobs 267271–267278 now publish shard collections,
build attachment in 18 independent tasks, index and audit. Historical details
below are retained for provenance.

The user authorized full distillation reconstruction on September 11, starting
with CIF databases. Earlier comments deferring full builds are superseded.

Controller job **266702** (resuming 266684) builds nine stages using the native DataCooker planner:
RNA, short, disordered and long base/attached CIF databases, including the final
disordered contact graph. Up to 18 worker nodes run each stage; the controller
uses node01 and worker submissions exclude node01 and node04. Node04 is left out
of this workload for experiments, per the updated user preference. This does not
create a cluster-wide reservation or displace other users. Merge, index and audit gates run after each stage. All heavy input
scans, sample recipe execution and audits take place on compute nodes.

Inputs are the provided `BioMol/materials/openfold_distillation` snapshot, the
existing clean CCD, sequence-ID map and 40% clustering. This run does not download
a different snapshot or run MSA/template processing. Outputs are under
`BioMol_clean/lmdb/distillation_{short,long,rna,disordered}/cif/`.

Source/configuration snapshots and their SHA256 manifest are retained in
`logs/distillation_cif/20260911/release/`. A 16 TiB sparse LMDB address-space ceiling
is used consistently for workers and merge, avoiding the default 300 GB merge
limit; this does not preallocate 16 TiB of physical storage or RAM.

Each stage checks input-key uniqueness, full output/index key coverage, metadata
counts and sizes, and decodes every payload to check structure and attachment
invariants. Base input samples pass a recipe gate before bulk submission. Final
payload samples include 128 seeded random records and the 16 largest entries.
Base samples are compared to source coordinates and atom names. Attached samples
are checked against supplied FASTA, sequence IDs and cluster references, and
against parent coordinates after water removal. Disordered contact samples are
checked independently with KD-tree neighbor counts. Available production attached
databases are compared for complete key coverage and sampled payload differences.
Differences are retained for analysis rather than silently called equivalent.

Report job **266685** runs after any controller exit, records its terminal state,
and writes `logs/distillation_cif/20260911/SUMMARY.md` and `completion.json`.
Failures stop subsequent stages and preserve source lists, missing-key reports,
shards and logs. Resuming accepts existing outputs only with matching run state and
a retained successful terminal-job record. Completed audits are reused; interrupted
audits restart without rebuilding the DB. Untracked outputs cause a stop.

The resource expansion replaced only the controller while it was auditing RNA;
both RNA builds were already complete and no worker was cancelled. Resource policy
is read from `logs/distillation_cif/20260911/resources.json` before each submission.

This is a submitted run, not a completion claim. No interactive polling loop is
used. Inspect the final report and production differences after execution before
declaring end-to-end completion.

## Recovery after audit dependency failure

Controller 266702 stopped during disordered contact auditing because the original
reference checker imported SciPy, which is absent from the project environment.
This was a checker dependency failure after DB generation, not a demonstrated
contact-data mismatch. Logs and the failed run reports are preserved.

Controller **266819** resumes these outputs with frozen orchestration code in
`logs/distillation_cif/20260911/scheduler_v2`. Its contact reference uses an independent
NumPy x-axis sweep. Nine regression tests passed on compute job **266818**.
Audits now run separately on node01 while the next native build progresses with
up to 18 workers. Report job **266821** follows controller termination. All audit
jobs must succeed before final completion; this recovery is not yet a completion claim.

## Direct SLURM execution after inventory bottleneck

The long dataset contains 16,099,484 entry directories. Recursive per-directory
inspection blocked the prior controller before worker submission. Controller
266819 was cancelled after the user explicitly authorized SLURM execution.
Inventory job 266822 enumerated directory names once and wrote 18 manifests;
NPZ existence/readability is checked by native build workers. Missing NPZs remain
coverage failures for investigation rather than silently disappearing from inputs.

Base array 266823 is released (nine tasks initially running; remaining tasks were
limited by QOSMaxCpuPerUserLimit). Merge/index: 266824 / 266825. Direct attached
array: 266848 (18 tasks, 56 CPUs each), merge/index: 266849 / 266850. Independent
base/attached audits: 266851 / 266852. Final report: 266853. No resident controller
is used by this chain. Pending old controller/report jobs 266826 / 266827 were
cancelled. The experiment exclusion is now node20; running base workers were
preserved and all still-pending base tasks explicitly exclude node20. Later
submissions embed exclusions in SLURM scripts/arguments rather than relying on
an ineffective SBATCH_EXCLUDE environment variable.

A separate disordered contact audit discrepancy was found after replacing SciPy.
It remains unresolved and is not treated as a successful audit; it does not block
the independent long build. The final long report does not claim completion of
all distillation CIF datasets. Job IDs and scripts are retained under
`logs/distillation_cif/20260911/long_run/`.
