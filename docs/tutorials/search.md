# Search (clustering · MSA · templates)

The search workflows turn extracted FASTAs into clusters, MSAs, and template
hits. They produce files, not LMDBs, so `structcooker build` runs them through
`datacooker.cli.workflow` as `materialize` (one job) or `parallel` (a split
recipe fans the items out; `node_count` makes it an N-node array) ops.

## Sequence clustering

mmseqs2 clustering at 30 % identity (with SAbDab antibody handling) → a cluster
table used by metadata attach, filtering, and template search. Clustering is
**per DB** while the sequence id space is shared (see the
[seq-id and cluster scheme](../seq-id-and-cluster-scheme.md)); every cluster id
is written as `c<DB>_<rep seq_id>`.

- Recipe: `structcooker.workflows.search.seq_cluster` (`cluster_dict`)
- Writer: `...writers.projections.write_clusters`
- Configs (materialize): `db/metadata/pdb_seq_cluster30.yaml` (`db_code: PDB`
  → `metadata/seq_cluster30_PDB.tsv`), `ofd_seq_cluster30.yaml` (`OFD`),
  `teddymer_seq_cluster30.yaml` (`TDM`), `afm_seq_cluster30.yaml` (`AFM`)
- Parameters (all four): `mmseqs2_seq_id: 0.30`, `mmseqs2_cov: 0.8`,
  `mmseqs2_covmode: "0"`, `mmseqs2_clustermode: "1"`; inputs are the DB's own
  FASTA, `seq_id_map.tsv`, and the SAbDab summary
  (`BIOMOL_ROOT/materials/raw/sabdab/sabdab_summary_all.tsv`).

```bash
pixi run structcooker build metadata/pdb_seq_cluster30
```

## MSA search

The PDB and OpenFold MSAs are provided inputs; the release build does not
rerun their HHblits search (the `search.msa_search_hhblits` recipe has no `db/`
config). See [MSA reproduction](../msa-reproduction.md).

Teddymer MSAs were searched once with **MMseqs2** (CPU) against UniRef30 2302 only; the
result is a supplied input (`materials/reference/msa_teddymer`) and the release build does
not search again. `teddymer/msa_search` is outside the release DAG, for sequences that have
no MSA yet:

- Split recipe: `structcooker.workflows.metadata.load_protein_sequence_chunks`
- Recipe: `structcooker.workflows.search.msa_search_mmseqs`
- Config (parallel): `db/teddymer/msa_search.yaml` — `threads: 28`,
  `chunk_size: 1000`, `n_jobs: 4` (4 × 28 = 112 CPUs per node), `node_count: 8`;
  queries from `materials/intermediate/worklists/seq_id_map_TDM.tsv`
- Output: `<seq_id>.a3m` under `BIOMOL_ROOT/materials/reference/msa_teddymer/`, in the
  layout the a3m ingest reads (`db/teddymer/msa.yaml`).

```bash
pixi run structcooker build teddymer/msa_search   # an 8-node array
```

## Template search

hmmbuild turns each query a3m into an HMM, and hmmsearch runs it against the
PDB polypeptide(L) chains. The earlier HHsuite path (hhmake + hhsearch) was
replaced by hmmsearch and removed.

| Stage | Recipe | Config |
| --- | --- | --- |
| target DB | `exports.pdb_polypeptide_l` | `db/metadata/pdb_polypeptide_L.yaml` (materialize → `pdb_polypeptide_L.fasta`) |
| strip insertions | `filters.a3m_cleanup` | `db/template/msa_wo_lower.yaml`, `db/teddymer/msa_wo_lower.yaml`, `db/afdb_multimer/msa_wo_lower.yaml` (parallel) |
| hmmbuild + hmmsearch | `search.hmmsearch` (`hmmsearch_results`, `hmmbuild_results`) | `db/template/hmmsearch.yaml` (PDB), `db/teddymer/hmmsearch.yaml`, `db/afdb_multimer/hmmsearch.yaml` (parallel, `node_count: 8`) |

hmmbuild needs a fixed-width alignment, so the queries are the
lowercase-stripped a3m files. The searches are idempotent (an existing HMM /
`.out` is skipped), and the teddymer / AFM configs list only pending queries
(`search.load_pending_a3m_list`) so a resume spreads the remainder evenly.

```bash
pixi run structcooker build template/hmmsearch    # node_count: 8 -> an 8-node array
pixi run structcooker build teddymer/hmmsearch
```

The hmmsearch output feeds [Template ingest](ingest-template.md).
