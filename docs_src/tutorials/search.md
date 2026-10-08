# Search (clustering · MSA · templates)

The search workflows turn the extracted FASTA into clusters, MSAs, and template
hits. They are compute-heavy and run through `datacooker.cli.workflow`
(`run` / `parallel-run`) rather than the LMDB builder.

## Sequence clustering

mmseqs2 clustering at 30 % identity (with SAbDab antibody handling) → a cluster
table used by metadata, filtering, and template search.

- Recipe: `structcooker.workflows.search.seq_cluster` (`cluster_dict`)
- Config: `configs/legacy/search/seq_cluster30.yaml` (`mmseqs2_seq_id: 0.3`, `cov: 0.8`)

```bash
pixi run python -m datacooker.cli.workflow run configs/legacy/search/seq_cluster30.yaml
```

## MSA search

hhblits against **UniRef30** and **BFD**, parallelised over sequence splits.

- Split recipe: `structcooker.workflows.metadata.load_protein_sequences`
- Recipe: `structcooker.workflows.search.msa_search` (`msa_results`)
- Config: `configs/legacy/search/msa_search_hhblits.yaml` (`db_ur30`, `db_bfd`, `cpu_per_job`, …)

```bash
sbatch submits/search/msa_search.sh
# directly (parallel across splits):
pixi run python -m datacooker.cli.workflow parallel-run configs/legacy/search/msa_search_hhblits.yaml
```

## Template search

Build profiles and search for templates over the MSAs. Several stages:

| Stage | Recipe | Config |
| --- | --- | --- |
| hmmbuild + hmmsearch | `search.hmmsearch` (`hmmsearch_results`, `hmmbuild_results`) | `db/template/hmmsearch.yaml` (PDB), `db/teddymer/hmmsearch.yaml`, `db/afdb_multimer/hmmsearch.yaml` |

The earlier HHsuite path (hhmake + hhsearch) was replaced by hmmsearch and removed.

```bash
pixi run structcooker build teddymer/hmmsearch   # node_count: 8 -> an 8-node array
```

The hmmsearch output feeds [Template ingest](ingest-template.md).
