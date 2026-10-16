# Licenses

The databases here are derived from public sources, each under its own terms. Anyone who
redistributes this tree, or a database cut from it, carries those terms with it: credit the
sources below and keep this page with the data. The code that built it (StructCooker) is
MIT-licensed; that licence covers the code, not the data.

**Teddymer is non-commercial (CC BY-NC 4.0).** `lmdb/teddymer/*` and the tables derived from it
carry that restriction, so a copy of this tree that includes them is not freely usable for
commercial work. Everything else below permits commercial use with attribution.

Terms were read from each source's own page on 2026-10-16.

| source | used in | terms (as published by the source) | where |
|---|---|---|---|
| PDB (wwPDB mmCIF) | `lmdb/pdb/*`, `lmdb/chain/*`, PDB MSAs and templates | CC0 1.0: no restriction; citing the PDB entries and wwPDB is requested | https://www.wwpdb.org/about/usage-policies |
| PDB Chemical Component Dictionary | `lmdb/ccd/*` | CC0 1.0 (part of the wwPDB archive) | https://www.wwpdb.org/data/ccd |
| OpenFold3 training data (distillation sets) | `lmdb/distillation_*/*` | **CC BY 4.0**, one licence for the whole `openfold3-data` bucket; no further restriction stated. Requested citation: "OpenFold3 Training Data was accessed on DATE from https://registry.opendata.aws/openfold3", plus the OpenFold3 repository (https://doi.org/10.5281/zenodo.19001000) and manuscript (https://www.nature.com/articles/s41592-024-02272-z). The sets are built from PDB and MGnify; the page does not restate their terms. | https://registry.opendata.aws/openfold3 |
| AlphaFold DB (multimer set) | `lmdb/afdb_multimer/*` | **CC BY 4.0** for AlphaFold DB (academic and commercial use, attribution; cite Jumper et al., Nature 2021). The complex release is served from EBI's `collaborations/nvda/` area, and EBI lists separate notices for some collaboration datasets: **confirm the notice for this release** before publishing. | https://alphafold.ebi.ac.uk/assets/License-Disclaimer.pdf , https://alphafold.ebi.ac.uk/download |
| Teddymer | `lmdb/teddymer/*`, `metadata/seq_cluster30_TDM.tsv` | **CC BY-NC 4.0**: "All files are available under a Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)". Non-commercial use only; attribution required, no citation text is given on the page. Whether the restriction extends to models trained on it is not settled here. | https://teddymer.steineggerlab.workers.dev/ |
| SAbDab | antibody chain split in the clustering | **CC BY 4.0** ("the data contained within SAbDab is available under a CC-BY 4.0 license", OPIG). Read from the OPIG pages as search results quote them: the current site renders in JavaScript, so re-check in a browser. | https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab |
| UniRef30 (via the Teddymer MSAs) | `lmdb/teddymer/msa*` | UniProt data: CC BY 4.0 | https://www.uniprot.org/help/license |

## Where the Teddymer terms reach

Files built from Teddymer: `lmdb/teddymer/*` (structures, MSAs, templates) and
`metadata/seq_cluster30_TDM.tsv`. Entries from Teddymer are also present in the tables shared
across sources, `metadata/seq_id_map.tsv` and the per-source and global cluster tables
(`seq_cluster30.tsv`, `seq_cluster40.tsv`), and in `materials/{raw,reference}/teddymer*`.
A copy meant for unrestricted or commercial use leaves out `lmdb/teddymer/*`, the TDM table and
the Teddymer rows of the shared tables. Those rows are the ids listed in
`materials/intermediate/worklists/seq_id_map_TDM.tsv` (999,853 ids; ids carry no source prefix).

`docs/scenarios.md` lists which uses need Teddymer; the PDB, PDB + distillation and FoldBench
scenarios do not.
