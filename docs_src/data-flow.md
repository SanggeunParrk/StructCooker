# Data flow

The StructCooker workflows compose into one end-to-end pipeline that turns raw
wwPDB / sequence data into the LMDB databases used downstream. Each box is a
workflow with its own [tutorial](tutorials.md); arrows are the artifacts they
produce and consume.

```mermaid
flowchart TB
    ccd[CCD *.cif] -->|ingest/ccd| ccddb[(CCD LMDB)]
    cif[RCSB *.cif] -->|ingest/cif + CCD LMDB| cifdb[(CIF LMDB)]
    cifdb -->|exports/fasta| fasta[merged.fasta]
    fasta -->|metadata/seq_id_map| sid[seq_id_map.tsv]
    fasta -->|search/seq_cluster · mmseqs2| clu[seq_cluster30.tsv]
    fasta -->|search/msa_search · hhblits| a3m[a3m files]
    a3m -->|search/hmmsearch + hhsearch| hmm[template hits]
    a3m -->|ingest/a3m| msadb[(MSA LMDB)]
    hmm -->|ingest/template_lmdb + CIF LMDB| tmpldb[(template LMDB)]

    cifdb -.->|filters/data · train/valid| split[(train / valid LMDB)]
    split -.->|metadata/attach| att[(attached LMDB)]
    cifdb -.->|exports/interacting_seq_*| inter[interacting clusters]
```

## Main flow

1. **CCD LMDB** — ingest the wwPDB Chemical Component Dictionary
   ([CCD ingest](tutorials/ingest-ccd.md)).
2. **CIF LMDB** — ingest RCSB mmCIF structures, resolving components against the
   CCD LMDB ([CIF ingest](tutorials/ingest-cif.md)).
3. **FASTA** — extract sequences from the CIF LMDB
   ([Exports](tutorials/exports.md)).
4. **Sequence id map** — assign stable sequence ids
   ([Metadata](tutorials/metadata.md)).
5. **Sequence clustering** — mmseqs2 clusters at 30 % identity
   ([Search](tutorials/search.md)).
6. **MSA search** — hhblits against UniRef30 + BFD
   ([Search](tutorials/search.md)).
7. **Template search** — hmmsearch / hhsearch over the MSAs
   ([Search](tutorials/search.md)).
8. **MSA LMDB** — pack the alignments, keyed by sequence id
   ([MSA ingest](tutorials/ingest-msa.md)).
9. **Template LMDB** — resolve template hits into `TemplateMol`s against the CIF
   LMDB, keyed by sequence id ([Template ingest](tutorials/ingest-template.md)).

## Optional flow

- **Filter** the dataset into train / valid splits
  ([Filters](tutorials/filters.md)).
- **Attach metadata** to the filtered structures
  ([Metadata](tutorials/metadata.md)).
- **Extract** interacting sequence ids / clusters and other artifacts
  ([Exports](tutorials/exports.md)).
- **Profile / validate** the resulting databases
  ([Analysis](tutorials/analysis.md)).

!!! note "Engine vs. domain"
    Every box above is a [DataCooker](https://CSSB-SNU.github.io/DataCooker/)
    recipe behind a config. The *engine* mechanics (recipes, `cli.lmdb`,
    `cli.workflow`) are documented in the DataCooker docs; these tutorials cover
    the *domain* recipes and how to run them.
