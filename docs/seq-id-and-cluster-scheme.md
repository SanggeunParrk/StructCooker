# Sequence identity and clustering — one id space, per-DB clusters

**Decided 2026-09-23.** Supersedes the single global clustering that produced
`metadata/seq_cluster30.tsv` / `seq_cluster40.tsv`.

## The rule

| | Scope | Why |
|---|---|---|
| `seq_id` | **one space, shared by every DB** | The id names a *sequence*. Two DBs holding the same sequence must agree on its id or nothing can be joined across them. |
| `seq_cluster` | **one clustering per DB** | The cluster names a *neighbourhood*, and what counts as a neighbourhood depends on what the DB contains. |

A cluster id carries the DB it was computed in:

    c{DB}_{seq_id}          e.g. cPDB_P00000000000000512758

## DB codes

Three letters, fixed width, so the code is a slice and ids line up in columns.

| Code | DB | State |
|---|---|---|
| `PDB` | PDB | base + attached built |
| `OFD` | OpenFold distillation (short / long / rna / disordered) | base + attached built |
| `TDM` | Teddymer | base built |
| `AFM` | AFDB multimer | raw only — no LMDB yet |

`AFD` (AFDB monomer) is deliberately **not** reserved: the raw tarballs are on disk but
there is no recipe and no plan. Adding a code later touches nothing that already exists.

`TDM`, not `TED`, for Teddymer: TED (The Encyclopedia of Domains) is Teddymer's upstream
and ships as its own database (`teddb.tar.gz`). Keeping `TED` free avoids the collision
if it is ever ingested.

## Why clustering is per-DB

One global clustering is what production did, and it makes clusters that do not mean
anything in particular. Measured on the existing `seq_cluster40.tsv`:

- 4,819,226 clusters over 16,902,243 members, **78.3% of them singletons**
- distillation sequences are in there, nearly all as singletons (`cP...` = their own id)

That is the predictable outcome of clustering unrelated populations together. PDB chains
are experimentally determined, mostly from a few well-studied organisms; OpenFold
distillation sequences are MGnify metagenomic predictions; Teddymer chains are not whole
proteins at all but TED domains cut out of AFDB models. A 40%-identity neighbourhood
computed across all three answers no question anyone asks: a PDB chain's nearest
metagenomic neighbour does not make it redundant for PDB-based splits, and a domain does
not cluster with the whole proteins it was cut from.

Per-DB clustering makes each cluster answer the question its DB is used for -- "what else
in *this* DB is redundant with this sequence" -- which is what train/valid splitting and
sampling weights actually need.

The shared `seq_id` keeps cross-DB joins available: to ask whether a Teddymer domain
sequence also appears in PDB, compare seq_ids, not clusters.

## What this changes

Cluster-id construction, six places:

| Where | Now |
|---|---|
| `instructions/transforms/metadata.py:164`, `:331` | `f"c{seqid}"` when a seq_id is absent from the map |
| `instructions/transforms/graph.py:812-813` | same pattern, interface clusters |
| `instructions/writers/projections.py:36`, `:82` | `f"c{rep_seq_hash}"` when writing the cluster TSV |

All of them need the DB code as an input, which means the recipes must carry one. The
functions currently have no way to know which DB they are processing.

Unaffected: `utils/mapping.py`'s `cluster_types` / `cluster_maps` map the *molecule type*
letter (`P`/`R`/`D`/...). Under the new id that letter is the first character after `_`.

`seq_id` itself does not change. Its rule stays as
`instructions/transforms/sequence.py:build_seq_id_map` defines it: a running counter per
molecule type, seeded from a provided map so existing ids survive a rebuild.

## Per-DB FASTA

Each DB gets exactly one consolidated FASTA, which is both the clustering input and the
attach input:

| DB | FASTA | Source |
|---|---|---|
| `PDB` | `materials/raw/fasta/cif_pdb.fasta` | extracted from `pdb/cif` |
| `OFD` | `materials/raw/fasta/ofd.fasta` | **to build** — currently four separate provided files (short / long / rna / disordered) |
| `TDM` | `materials/raw/fasta/teddymer.fasta` | extracted from `teddymer/cif` (recipe exists, not yet run) |

Splitting OFD across four files is why `seq_id_map.yaml` lists four paths; one file per DB
makes the clustering input and the DB a one-to-one thing.

## Migration

Every attached DB currently holds old-format ids (`cP...`), PDB and distillation alike, so
the format alone cannot say which DB an id came from. There is no compatibility rule to
write; the attached DBs are re-attached.

Re-attaching does **not** require rebuilding the base CIF LMDBs -- attach reads a base and
writes a new attached DB. `BioMol_clean` holds every base (production does not; it kept
only the attached ones). Measured read+write:

| DB | Entries | IO |
|---|---:|---:|
| distillation_long | 16,099,484 | ~2,085 GB |
| PDB | 249,676 | ~233 GB |
| teddymer | 510,454 | ~43 GB |
| distillation_short | 430,418 | ~24 GB |
| distillation_rna | 126,778 | ~8 GB |
| distillation_disordered | 28,567 | ~6 GB |

Long is 87% of the work; everything else together is a few hours. Attach is IO-bound (it
adds two chain fields), so cores do not help -- disk bandwidth sets the pace.

## Open

- `reference/seqcluster_corpus.fasta` points at a **PDB-only** 6,219,105-sequence file
  (`BioMolDB_20260224/fasta/merged.fasta`) while the production cluster file covers all
  16,902,243 seq_ids. The union corpus production actually clustered is
  `openfold_distillation/fasta/all.fasta` -- measured 2026-09-23: 23,015,909 sequences =
  the four distillation files (16,796,804 together) plus exactly those same 6,219,105 PDB
  entries. The reference pointer holds half the corpus. Per-DB clustering removes the need
  for a union corpus at all, so retire the pointer rather than repoint it.
- `pdb/cif_pdb_attached` holds 32 fewer entries than `pdb/cif` (249,644 vs 249,676). The
  re-attach should surface why.
- `distillation_long` has stale `.build.lock` files and an
  `.cif_long_attached.empty_attempt_267274.lmdb` (0 entries) from a failed attach. Check
  before re-attaching.
