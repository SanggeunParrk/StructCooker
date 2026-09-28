# Long CIF recovery — September 14

The original 18 base shards and 18 attached shards each contained 16,099,470
records. Native attached reports recorded zero failed records. Both coverage
audits correctly failed against the 16,099,484-entry input inventory.

All 14 omissions failed with `Bad magic number for central directory`. Each NPZ
contained a complete, valid ZIP archive followed by 32–479 stale trailing bytes,
including a later end-of-directory marker. The parser selected that later marker.
Recovery retained the valid archive prefix in separate files; shared originals
were never overwritten. All 13 ZIP members passed CRC checks in every recovered
file, and all 14 files passed the actual structure-conversion recipe.

Evidence is retained in `logs/distillation_cif/20260914/`:

- `missing_npz_diagnostics.json`: original hashes, ZIP errors and header offsets.
- `npz_recovery.json`: original/recovered hashes, retained lengths, member names
  and successful recipe validation.
- `recovered/`: exact valid-prefix NPZ copies used to build the additional shard.
- `patch_build_complete.json`: 14 base and 14 attached records built successfully
  and structurally checked, with native failure ledgers alongside both databases.

The same 14 entries were absent from the supplied FASTA and sequence-ID reference.
Their FASTA sequences were extracted from the recovered structures. None matched
the supplied sequence-ID map. A separate reference version preserves all existing
IDs and appends 14 protein IDs above the existing maximum, assigned deterministically
by sequence order. Existing reference files remain unchanged:

`BioMol_clean/metadata/long_recovered_20260914/{seq_id_map.tsv,long.fasta}`

`metadata_version.json` records the exact new IDs and reference paths. The new IDs
are absent from the supplied clustering, so attachment uses the existing recipe's
singleton fallback `c<seq_id>`. This is not a claim that the recovered sequences
were reclustered at 40% identity. Consumers of the recovered release, including
later MSA key generation, must use this reference version instead of the older
reference files. Do not combine ID namespaces from independently extended versions.

## Publication and verification

No bulk payload merge or rebuild is required. Each logical DB receives one additional
14-record shard, bringing it to **19 shards and 16,099,484 records**. The additional
shards are copied into `.long_recovered_20260914/` under the clean long CIF directory.
Publication checks that their keys are disjoint from the existing collection,
preserves old manifests/indexes, and extends indexes using existing size rows.

| Step | SLURM job |
|---|---|
| Inspect damaged inputs | 267339 |
| CRC and conversion verification | 267340 |
| Completed base/attached patch and reference version | 267347 |
| Reconcile attached publication (base patch already published) | 267375 |
| Full paired shard audits on node01–19 | 267353–267371 |
| Aggregate attached index and report | 267376 |
| Small audit of all 14 repaired records | 267379 |

Each audit checks exact expected key coverage in its base and attached shard,
index sizes and every decoded structure's invariants. Seeded and largest-record
samples are compared to original NPZ coordinates/atom names, parent coordinates
after water removal, and the versioned sequence/cluster references. All 14 recovered
records are included in these source and metadata comparisons. Production sample
differences are recorded separately; the recovered records may be absent from
production. Node20 remains excluded.

Audit source code is frozen under `audit_release/`. Individual `audit_<i>.json`
files retain progress. The final `completion.json` distinguishes failed checks
from `audited_pending_analysis`; it does not establish MSA/template completion or
resolve the independent disordered-contact discrepancy. Inspect completed reports
after the jobs finish instead of continuously polling.

## Publication reconciliation

The earlier queued attached publisher retained the initial, empty shard-list path
in SLURM's submitted script, despite the local script later being edited. Job
267274 consequently published an empty LMDB. The completed `.long_attached_shard_policy_v2`
shards were intact; their final ledgers and actual entry counts sum to 16,099,470.
This was a publication error, not a loss of those records.

Job 267352 successfully extended the base manifest/index, then stopped when it
encountered that empty attached LMDB. Job 267375 verifies the correct 18 attached
shards plus the recovery shard and preserves the empty artifact separately before
publishing the correct manifest. The audit jobs were reconnected to 267375. Each
audit also creates its attached physical-shard index; final job 267376 consolidates
those indexes only if all audits pass. Superseded report job 267372 was cancelled.

## Completed audit analysis — September 15

Both databases passed all 19 shard audits, covering 16,099,484 records each.
The production sample cluster differences have been traced to different 40%
reference snapshots; 29 of 32 compared cluster member sets changed and three
retained their members under different labels. See
[the follow-up report](distillation-followup-2026-09-15.md) for evidence and
the remaining contacts/MSA/template jobs.

## ID namespace resolution — September 26

The per-DB seq_id/cluster rebuild ([scheme](seq-id-and-cluster-scheme.md)) regenerated
`BioMol_clean/metadata/seq_id_map.tsv` from the production reference, independently of the
recovery version above. Its re-attach of long therefore used the supplied `long.fasta`,
which lacks the 14 recovered entries, and published 16,099,470 records. The recovery-version
IDs `P…16902258`–`P…16902271` could not be reused: the regenerated map had already issued
those numbers to other sequences — the collision this document warned against.

Resolution: the recovery version is superseded for seq_id purposes.

- The 14 records (the bytes after the supplied file's end in the recovered `long.fasta`;
  prefix verified byte-identical) are `BioMol_clean/materials/raw/fasta/ofd_long_recovered.fasta`,
  appended **last** in `db/metadata/seq_id_map.yaml` so no earlier ID moves.
- They received `P00000000000019613388`–`P00000000000019613401`. Verified (job 274028):
  19,613,276 existing IDs kept, 0 changed, 0 lost, 14 new, all 14 equal to the recovered sequences.
  Pre-change map: `metadata/seq_id_map.before_long14.tsv`.
- `db/distillation/long_cif_attached.yaml` reads the recovered `long.fasta` (supplied + 14).
  They are not in `seq_cluster30_OFD.tsv`, so they take the singleton fallback `cOFD_<seq_id>`,
  as in the September 14 recovery. The 16,099,470-record attached DB is kept in
  `lmdb/distillation_long/cif/_superseded_20260926/`; rebuild jobs 274029–274031.
