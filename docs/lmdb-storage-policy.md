# LMDB storage policy — 2026-09-13

September 14: the 14 omitted long CIF entries were recovered from CRC-verified NPZ
prefixes. The [recovery report](long-cif-recovery-2026-09-14.md) describes the extra
shards, required metadata reference version, and queued full audits. The table
below records the pre-recovery inventory; the recovered long CIF target count is
16,099,484 per DB.

After all build shards finish, count their actual output records. At **5,000,000
items or fewer**, merge into one LMDB. Above **5,000,000**, publish a `shard_lmdb`
manifest and retain the physical shards as permanent data. MSA depth is not the
item count. A chain-exploding recipe uses its output chain count, not input CIF count.

The finalizer checks the threshold at runtime. The table is an inventory, not a
hard-coded dataset allowlist. Existing merged databases are not split or rebuilt
just because the policy changed.

| Dataset | Variant | Items | Storage |
|---|---|---:|---|
| CCD | CCD | 50,782 | merged |
| PDB CIF | base / attached | 249,676 / 249,644 | merged |
| PDB chain | CIF | 4,111,034 | merged |
| PDB train | 20210930 / 20260301 / MPNN | 167,912 / 233,579 / 484,383 | merged |
| PDB valid | valid1 / attached / valid2 | 21,527 / 21,527 / 2,878 | merged |
| PDB protein MSA | base / d16k / d2k / d512 | 178,249 each | merged |
| PDB RNA MSA | base / d2k / d512 | 6,572 each | merged |
| PDB template | seqid / chain | 176,579 / 1,007,928 (reference) | merged |
| Distillation long CIF | base / attached | 16,099,470 (base shards / attached reference) | shard_lmdb |
| Distillation long MSA | base | 6,711,214 already generated, partial | shard_lmdb |
| Distillation long MSA | d2k / d512 | 16,098,796 each (reference) | shard_lmdb |
| Distillation long template | topn | 16,080,759 (reference) | shard_lmdb |
| Distillation short CIF | base / attached | 430,418 each | merged |
| Distillation short MSA | base / depth caps | 430,245 (caps expected equal) | merged |
| Distillation short template | topn | 430,370 (reference) | merged |
| Distillation RNA CIF | base / attached | 126,778 each | merged |
| Distillation RNA MSA | base / depth caps | 126,751 (caps expected equal) | merged |
| Distillation disordered CIF | base / attached / contacts | 28,567 each | merged |
| Distillation disordered MSA | base / d2k / d512 | 19,649 (caps reference; base expected equal) | merged |
| Distillation disordered template | template | 63,610 (reference) | merged |

Counts were read from LMDB metadata on September 13. Reference rows are existing
production outputs, not claims that new builds have completed. PDB chain template
uses the existing `template_reldate.lmdb` reference. TSV/FASTA projections are outside
this storage rule.

## Format and readers

The logical target remains a directory named `<name>.lmdb`. A sharded target has
`shards.json` instead of `data.mdb`. The versioned manifest stores relative physical
shard paths and counts. Publication scans keys to reject duplicates, then atomically
replaces the manifest; it never copies payloads. A large build with colliding keys
fails instead of silently applying overwrite/union semantics across shards.

Use `datacooker.lmdb.sharded.open_env(path, readonly=True, lock=False)` for either
format. DataCooker's count, extraction, rebuild and indexing paths use this reader.
It supports get, entry count and globally sorted forward iteration. It is not a
replacement for every python-lmdb API: direct third-party `lmdb.open` callers need
to adopt the reader. Random lookup currently probes shards; partition-aligned
rebuild jobs should open their corresponding physical source shard directly.

The logical DB retains `.index.tsv` and `.meta.json` sidecars. Indexing scans up to
eight physical shards concurrently. Content fingerprints include the manifest and
every referenced shard; manifest existence alone does not establish verified completion.

Treat referenced shard directories as permanent: do not clean up a `.build` or
stream-build directory while a manifest references it. Keep shards immutable and
preserve their relative layout when relocating a collection. Never replace a live
physical shard in place. This change avoids the merge copy and its temporary disk
duplication, but does not eliminate the work of indexing and auditing.

## Long CIF migration

Replacement scripts, frozen source, and job IDs are under
`logs/distillation_cif/20260913/shard_policy/`. The 18 completed base shards are
preserved. Attached workers read each physical base shard directly, so base
indexing/auditing and attachment can progress independently. Both outputs publish
manifests; neither schedules a payload merge. Node20 remains excluded.

The old partial merged output is not adopted or deleted by this migration. Its
staging directory remains available for an explicit cleanup after replacement
validation. The 14-entry discrepancy between long input inventory and built base
shards remains subject to coverage auditing; publishing a manifest does not waive it.

Replacement SLURM jobs submitted September 13:

| Stage | Job |
|---|---|
| Base shard publication | 267271 |
| Base index | 267272 |
| Attached build, 18 tasks | 267273 |
| Attached shard publication | 267274 |
| Attached index | 267275 |
| Base audit | 267276 |
| Attached audit | 267277 |
| Final report | 267278 |

Retired jobs: 266824, 266825, 266848, 266849, 266850, 266851, 266852, 266853.
No completed base-worker output was deleted. Full regression job 267270 passed
173 tests, including native rebuild from a shard collection. Ruff and Pyright
passed for both source packages. These are software checks, not a completed
16-million-record dataset audit.

Initial attachment startup exposed nested `Path` values being serialized as plain
strings. The worker now preserves the original OmegaConf path resolvers, and the
official production serializer has a regression-tested fix for the same issue.
The 18-task array 267273 was requeued into fresh `.long_attached_shard_policy_v2`
outputs; initial failure logs are archived under `initial_path_failure/`. The final
regression job 267296 passed 174 tests. `migration.json` records the replacement and
retry details. This retry does not regenerate base structures.
