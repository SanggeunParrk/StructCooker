# Tutorials

These tutorials cover the **domain workflows**. For how the recipe engine and
LMDB builds work, see the
[DataCooker docs](https://CSSB-SNU.github.io/DataCooker/). The
[data flow](data-flow.md) page shows how they fit together end to end.

## Ingest (build LMDBs)

- [CCD ingest](tutorials/ingest-ccd.md) — Chemical Component Dictionary.
- [CIF ingest](tutorials/ingest-cif.md) — RCSB mmCIF structures → `CIFMol`.
- [MSA ingest](tutorials/ingest-msa.md) — a3m alignments.
- [Template ingest](tutorials/ingest-template.md) — template hits → `TemplateMol`.
- [OpenFold3 distillation ingest](tutorials/openfold-distillation.md) —
  MSA / structure / template from the distillation datasets.

## Build the inputs

- [Search](tutorials/search.md) — sequence clustering, MSA search, template search.
- [Metadata](tutorials/metadata.md) — sequence id map, attach metadata.

## Curate & inspect

- [Filters](tutorials/filters.md) — train / valid splits, MSA filtering.
- [Exports](tutorials/exports.md) — FASTA, interacting sequences, graph tables.
- [Analysis](tutorials/analysis.md) — profiling and validation.
