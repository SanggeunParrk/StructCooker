# Licenses

The databases here are derived from public sources, each under its own terms. Anyone who
redistributes this tree, or a database cut from it, carries those terms with it: credit the
sources below and keep this page with the data. The code that built it (StructCooker) is
MIT-licensed; that licence covers the code, not the data.

| source | used in | terms (as published by the source) | check |
|---|---|---|---|
| PDB (wwPDB mmCIF) | `lmdb/pdb/*`, `lmdb/chain/*`, PDB MSAs and templates | CC0 1.0: no restriction; citing the PDB entries and wwPDB is requested | https://www.wwpdb.org/about/usage-policies |
| PDB Chemical Component Dictionary | `lmdb/ccd/*` | CC0 1.0 (part of the wwPDB archive) | https://www.wwpdb.org/data/ccd |
| AlphaFold DB (multimer set) | `lmdb/afdb_multimer/*` | CC BY 4.0: attribution required | https://alphafold.ebi.ac.uk/ (terms of use) |
| OpenFold3 distillation sets | `lmdb/distillation_*/*` | to confirm on the dataset page before redistribution | https://registry.opendata.aws/openfold/ |
| Teddymer | `lmdb/teddymer/*` | to confirm with the source before redistribution | https://teddymer.steineggerlab.workers.dev/ |
| SAbDab | antibody chain split in the clustering | to confirm with the source before redistribution | https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab |
| UniRef30 (via the Teddymer MSAs) | `lmdb/teddymer/msa*` | UniProt data: CC BY 4.0 | https://www.uniprot.org/help/license |

Rows marked "to confirm" are not settled here: read the source's current terms before
publishing data derived from it.
