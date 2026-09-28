# MSA ingestion and verification

The PDB MSA workflow consumes provided A3M alignments. It does not regenerate
HHblits or RNA search results. The default input is
`DATA_ROOT/BioMol/materials/intermediate/msa`, overridable with `MSA_ROOT`:

```text
P/<last-three-digits>/<sequence-id>/<sequence-id>.a3m
R/<last-three-digits>/<sequence-id>/<sequence-id>.a3m
```

The structured file patterns select final alignments and avoid HHblits/SignalP
intermediate directories. All seven outputs live under `OUTPUT_ROOT/lmdb/pdb/msa/`:
`a3m.lmdb`, `a3m_d16k.lmdb`, `a3m_d2k.lmdb`, `a3m_d512.lmdb`, `msa_rna.lmdb`, `msa_rna_d2k.lmdb`, and `msa_rna_d512.lmdb`.
The protein and RNA caps retain the query and first hits, then recompute the profile
and deletion mean. They do not select a random subset.

From a compute allocation, use the native workflow:

```bash
pixi run structcooker inspect --manifest db/MANIFEST_msa.yaml --strict
export DATACOOKER_MAX_NODES=2
pixi run structcooker build-all --manifest db/MANIFEST_msa.yaml
```

`DATACOOKER_MAX_NODES` limits each pipeline's worker nodes; independent pipelines can
overlap, so it is not a global cluster quota. The September 10 run sequences the
protein builds and uses one RNA worker node, preserving experiment capacity.

## Input versions

The provided RNA A3Ms reproduce all sequence features in the retained
`msa_rna_pdb_shard0_orig_20260709.lmdb` snapshot in all **6,572 records**. They do
not reproduce all contents of the later `msa_rna_pdb_shard0.lmdb`: **3,604 records**
have different alignments or depths in that later DB. Both snapshots and the new
DB have exactly the same keys; 2,833 records have corrected header metadata
relative to the older snapshot. Some provided
A3Ms contain only the query while the later DB contains additional hits. Rebuilding
from provided files must not be described as reproducing those later search results.
The production DB is retained unchanged.

## Corrections

- Header parsing resets the match for every row and retains the complete fallback
  ID. Unknown hits no longer inherit the preceding hit's species and identifiers.
- A3M insertion counts follow aligned columns; trailing lowercase insertions have
  no following column and are excluded from the deletion matrix. Query insertions
  do not change the aligned width. Invalid widths fail with a specific row number.
- Profiles use bounded counting chunks, avoiding a depth × length × alphabet
  one-hot allocation. Output dtypes and the existing deletion-count clipping remain.
- Depth caps must retain at least one query row.
- The common codec can read legacy Zstandard frames without a recorded content
  size and rejects truncated frames or inconsistent array payload lengths.

Before the full build, all sequence fields of 48 protein samples (random plus large
records) matched production; 13 had corrected header metadata. Four OpenFold
categories (short, long, RNA, disordered), three real samples each, passed native
recipe and independent statistics checks, including 512/2000-row caps. Those are
recipe checks, not full OpenFold DB/key-space reproduction or a bulk rebuild.

All seven databases have been built, indexed, and audited. Protein controller job
266583, its final index job 266611, and final audit job 266572 completed with exit
code 0:0. Progress and audit artifacts are under `logs/msa/`.

| Output | Records | Maximum rows |
| --- | ---: | ---: |
| `a3m.lmdb` | 178,249 | Uncapped |
| `a3m_d16k.lmdb` | 178,249 | 16,000 |
| `a3m_d2k.lmdb` | 178,249 | 2,000 |
| `a3m_d512.lmdb` | 178,249 | 512 |
| `msa_rna.lmdb` | 6,572 | Uncapped |
| `msa_rna_d2k.lmdb` | 6,572 | 2,000 |
| `msa_rna_d512.lmdb` | 6,572 | 512 |

The full source test suite passed 148 tests on a compute node (job 266582).
Ruff and Pyright passed; Pyright reported zero errors and warnings in 169 files.

All 178,249 protein input keys are present, with no unexpected keys. All four
protein databases also exactly cover their corresponding production key sets.
Every index key and metadata count/byte summary was checked. In each database,
144 records (128 seeded random plus 16 largest) passed independent profile and
deletion-mean checks; capped records also matched their parent query, rows, and
headers. All sequence fields in these samples equal production. In each database,
34 sampled records have header-only corrections. This is sampled payload
verification, not a claim that every protein payload was compared. See
`logs/msa/final_audit.json` and `logs/msa/audit_terminal_state.txt`.

A 10,000 × 256 synthetic MSA benchmark (three runs per implementation, separate
processes on one compute node) produced identical feature bytes. Median parse time
was 1.635 s before and 1.520 s after; median process peak RSS was 371.7 MiB before
and 67.2 MiB after. This measures parsing, not full pipeline throughput.

RNA base, 2,000-row, and 512-row outputs each contain 6,572 records. Full key and
index coverage passed; both caps were checked for every record against their parent
and independent statistics. See `logs/msa/rna_final_audit.json` (job 266573).

## Coverage against the recovered CIF snapshot

The recovered CIF FASTA contains 178,611 distinct L-protein sequences and 6,584 RNA
sequences, all resolved in the sequence-ID map. Provided final alignments cover
178,247 of those protein sequences and 6,572 RNA sequences. Two additional supplied
protein alignment keys are outside this CIF sequence set.

The missing 364 protein sequences comprise 13 all-X queries and 351 queries with
known residues. All 12 missing RNA sequences contain only N. Thus the actionable
search gap is 351 protein sequences; unknown-only sequences contain no residues to
search. Some missing protein directories contain intermediate HHblits results but
no final alignment. Intermediate files are not silently promoted to final MSAs.

The missing IDs and exact queries are retained in `logs/msa/cif_P_missing_msa.txt`,
`logs/msa/cif_R_missing_msa.txt`, and `logs/msa/missing_protein_search.tsv`. Completing
these searches is distinct from reconstructing the supplied alignment databases.
The final audit verified the protein coverage key hash against the completed DB
keys; `logs/msa/cif_coverage.json` records `verified_built_keys`.
