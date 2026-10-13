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

| Code | DB | Cluster table |
|---|---|---|
| `PDB` | PDB | `metadata/seq_cluster30_PDB.tsv` |
| `OFD` | OpenFold distillation (short / long / rna / disordered) | `metadata/seq_cluster30_OFD.tsv` |
| `TDM` | Teddymer | `metadata/seq_cluster30_TDM.tsv` |
| `AFM` | AFDB multimer | `metadata/seq_cluster30_AFM.tsv` |

All four are clustered the same way: mmseqs2 at 30% identity, coverage 0.8
(`db/metadata/*_seq_cluster30.yaml`, see [stages](stages.md#sequence-ids-and-clustering)).

`AFD` (AFDB monomer) is deliberately **not** reserved: there is no recipe and no plan. Adding a code later touches nothing that already exists.

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

## seq_id

`seq_id` is a running counter per molecule type (`instructions/transforms/sequence.py:
build_seq_id_map`), seeded from `materials/reference/seq_id_map_seed.tsv` so existing ids
survive a rebuild; sequences not in the seed are appended in the order of
`db/metadata/seq_id_map.yaml`'s `fasta_paths`, which is why a new DB goes last. The molecule
type letter (`P`/`R`/`D`/...) is the first character of the id, and so the first character
after `_` in a cluster id.

## Per-DB FASTA

Each DB has exactly one consolidated FASTA, under `materials/intermediate/fasta/`, which is
both its clustering input and its attach input:

| DB | FASTA | Source |
|---|---|---|
| `PDB` | `cif_pdb.fasta` | extracted from `pdb/cif` |
| `OFD` | `ofd.fasta` | the four supplied reference fastas, concatenated (`metadata/ofd_fasta`) |
| `TDM` | `teddymer.fasta` | extracted from `teddymer/cif` |
| `AFM` | `afm.fasta` | extracted from both AFM cif DBs, concatenated (`metadata/afm_fasta`) |

Every attached DB was rebuilt with these cluster ids after the switch from the union
clustering (`seq_cluster30.tsv` / `seq_cluster40.tsv`); the old layout is described in
[data provenance](history/data-provenance.md).
