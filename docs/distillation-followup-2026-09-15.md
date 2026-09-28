# Distillation follow-up — September 15

This run continues the supplied distillation MSA and template builds. It does not
repeat the completed PDB build or the completed long CIF build. Heavy work is a
native SLURM dependency graph; no resident controller or interactive polling loop
is required. Node20 is excluded from every submitted job.

## Long CIF comparison resolved

The final long base and attached databases each contain 16,099,484 records in
19 physical shards. Their full structural/index checks passed. In the production
comparison samples, 150 records differed only in cluster IDs; the 14 recovered
records were absent from production. These are **sample comparison counts**, not
a claim that every production payload was compared byte-for-byte.

For 32 of the changed samples, all rebuilt assignments match
`BioMol_clean/metadata/seq_cluster40.tsv`, while all production assignments match
`BioMol/metadata/seq_cluster40.tsv`. A cross-snapshot membership comparison found
29 changed member sets and three identical member sets with different labels.
Thus this is a reference-version difference, sometimes including changed cluster
membership. Retain the rebuilt database's declared clean 40% reference. Do not
change its IDs merely to match the older production snapshot.

The 14 newly recovered sequences still use documented singleton fallback cluster
IDs; this run does not claim to have reclustered those sequences at 40% identity.
The separately recovered structures are also explicitly mapped into template
query-length loading with `OPENFOLD_STRUCTURE_OVERRIDES`. Shared raw files remain
untouched. Their recorded recovery SHA256 values passed the recipe gate.

## Contacts correction

Independent x-axis sweep calculations found that every nonempty disordered graph
in a 32-assembly sample counted each atom pair twice. Disordered attachment now
uses `chain_contacts_grid(..., count_atom_pairs_once=True)`. The default preserves
the existing PDB count convention, so completed PDB data is not rebuilt.

Job **267677** builds a fresh disordered contacts DB, checks every record against
the independent reference, and only then replaces the publication while preserving
its predecessor and sidecars. `contacts_complete.json` is the completion evidence;
a submitted or running job is not completion.

## Native job graph

Operational inputs, source snapshot, scripts and reports are under
`logs/distillation_followup/20260915/`. `jobs.json` contains the submitted IDs.

| Jobs | Work |
|---|---|
| 267674 | Unique sequence inventory, existing long-shard key coverage, missing input lists |
| 267678, 267680, 267681, 267710 | Seven real-input recipe gates, MSA statistics/cap checks, template shape/index checks, subprocess audit gate |
| 267682 array | Validate 14 existing long MSA shards and short/RNA MSAs; generate long caps; validate and reuse existing short/RNA caps |
| 267683 | Deduplicate disordered inputs by their intended keys |
| 267684 array | Build missing long MSAs and their depth caps |
| 267685 array | Build disordered MSAs and depth caps |
| 267686–267688 arrays | Long, short and disordered templates |
| 267689–267695 | Validate coverage, apply 5M publication policy, create indexes and metadata |
| 267696 | Aggregate terminal results, including missing stages and unavailable source counts |

Raw MSA depth is capped at 65,536; derived caps are 2,000 and 512, plus 5,000
for short/RNA. There are 18 partitions for new builds. Existing long partial
shards contain 6,711,214 records; they are admitted only after full payload
validation. Existing short/RNA caps are checked against the corresponding base
arrays, query and headers before reuse. Rebuilt caps receive the same comparison.
Input paths are deduplicated deterministically in the fixed inventory order by
sequence ID (MSA) or entry/chain ID (template).

New physical outputs have `.release_20260915` in their path. Existing long raw
shards and legacy cap files are retained. The publishers use actual output count:
**at most 5,000,000 records -> LMDB; more -> shard_lmdb**. A zero-size result,
duplicate keys, inconsistent accounting, a readable-source recipe failure, or
failed validation blocks publication. Unreadable source files receive a separate
ledger with the path, hash, size and reader error. Reports distinguish
`complete_with_source_exclusions` from complete input coverage. Such exclusions
still require analysis; they are not silently represented as successful records.

Five initial reuse tasks (7–11) could not create multiprocessing semaphores on a
node with exhausted shared memory. They were requeued after changing audits to
independent subprocesses; healthy running tasks were retained. The new audit
execution passed gate 267710. Pending workers were reduced to 32 CPUs to fit
currently available 32-CPU node slots; already-running 56-CPU allocations were
preserved. Scheduler overrides are recorded in `resource_updates.json`. No other
user's jobs were stopped.

## Validation and completion boundaries

The source changes passed Ruff, Pyright (zero errors/warnings), and all 177 tests
on compute nodes (267679 and 267713). The run uses the immutable source snapshot
recorded in `source.json`; the subsequent override-manifest type guard does not
change the recorded, valid manifest's behavior.

MSA audits decode all records and verify shapes, bounds, headers and recomputed
statistics. Cap audits additionally compare the retained arrays and query to the
base. Template audits verify record layout, the 20-hit cap for monomers, and
four-backbone-atom indexing. Empty template records are counted separately; the
underlying recipe may skip absent/unusable chain hits. These checks do not claim
full per-hit equivalence to the older production database.

Inspect `completion.json`, each `*_complete.json`, `receipts/`, `audits/`,
`exclusions/`, and failed-job logs after work terminates. Until then the remaining
MSA/template stages are in progress, not production-complete. A missing completion
report or `full_input_coverage: false` requires follow-up rather than a success claim.

## September 17 inspection

Live LMDB metadata showed 8,515,366 newly stored long MSA records plus the
6,711,214 reused records (15,226,580 / 16,099,468, approximately 94.6% stored).
Both depth caps had 12,447,376 records (approximately 77.3%). Long templates
had 5,022,478 / 16,099,484 records stored (approximately 31.2%). These are
storage counts, not terminal validated-completion percentages.

Six long MSA workers, all 18 long template workers, and one reuse/cap worker
remained running. Short and disordered template shard receipts were complete,
but their publishers failed: the run script incorrectly requested a fixed
100 TB LMDB address map. The map is now sized from the physical source data.
Replacement publisher jobs are 268373 and 268374. Disordered MSA validation
also incorrectly passed the template-only string sentinel as a numeric MSA
depth. This is fixed; 268375 resumes validation from existing build reports and
creates the caps without rebuilding the base, followed by publisher 268376.
The aggregate report dependency was updated to include these replacements.

Long MSA shard 3 remains blocked by one readable source whose recipe failed
(P00000000000005681270); it is not waived as an unreadable source. Short MSA
coverage is also 173 keys below the current inventory and must be reconciled
before its publisher can pass. These prevent an overall completion claim.

Node04 had 96 allocated CPUs but CPU load about 3, with approximately 326 GiB
used memory and 269 GiB shared memory; its active heavy Ray workers were about
three. Node12 similarly had a CPU load about 5, whereas node18 was about 96.
This is consistent with the existing memory-admission gate restricting
concurrency on some colocated jobs. Large data volume alone does not explain
the uneven throughput. No running build was interrupted during this inspection.
