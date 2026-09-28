# Raw materials — what the builds read, and where it came from

Two trees hold build inputs, and they are not the same kind of thing:

| tree | holds | written by |
|---|---|---|
| `DATA_ROOT/BioMol/materials/raw/` | downloaded or provided **source** data | staging scripts, once |
| `OUTPUT_ROOT/materials/` (`BioMol_clean`) | **derived** inputs: per-DB FASTAs, search intermediates | the pipeline |

Every staged source directory carries a `README.md` (what it is, why this subset) and a
`SOURCE.tsv` (where each file came from, its size). The layout below is the contract the
`db/*.yaml` configs rely on; a README in the directory itself is the detailed record.

## `BioMol/materials/raw/` — source data

```
raw/
├── ccd/                    components_20260803.cif.gz            wwPDB CCD snapshot
├── cif/                    <mid 2 chars>/<pdb id>.cif.gz         RCSB mmCIF mirror, 1,101 shards
├── fasta/                  cif.fasta · cif_pdb.fasta · pdb_polypeptide_L.fasta
├── rna_alignment_arrays/   <pdb>_<chain>.npz                     RNA alignment arrays
│
├── teddymer/                                                     TDM — 510,454 TED domain pairs
│   ├── README.md
│   ├── SOURCE.tsv          dimer_index · source path · bytes
│   └── pdb/<last 3 digits of DimerIndex>/<DimerIndex>DI_<acc>_<TED pair>.pdb
│
└── afdb_multimer/                                                AFM — AFDB complex release
    ├── README.md
    ├── homodimer/          1,750,755  high-confidence homodimers (manuscript set)
    │   ├── SOURCE.tsv      entity · original name · original path · bytes
    │   └── pdb/<last 3>/AF-<n>-model_v1.pdb
    ├── heterodimer/        80,248  passes_quality_threshold heterodimers
    │   ├── SOURCE.tsv
    │   ├── chain_msa.tsv   chain -> monomer MSA entity, content-verified
    │   └── cif/<last 3>/AF-<n>-model_v1.cif.gz
    ├── msa/                1,799,837  monomer MSAs the structures use
    │   ├── README.md
    │   ├── SOURCE.tsv      entity · release batch dir / tar member / generated · used by
    │   └── <last 3>/AF-<n>-msa_v1.a3m.zst
    └── msa_generated/      MSAs the release does not ship (1), with README
```

### teddymer

Selected, not the full release: 510,454 non-singleton cluster representatives that pass the
source paper's interface filter (length > 10, pAE < 10, pLDDT > 70) out of 10,089,503 dimers.
Each file is two TED domains of ONE AFDB v4 model — a domain-domain interface, not two
proteins. A domain can be sequence-discontinuous, so a chain's residue numbers may jump; the
jump is a domain boundary, not missing structure. B-factor holds pLDDT.
Staged by `scripts/maintenance/stage_teddymer.py` (copy).

### afdb_multimer

- **Homodimers** are the manuscript bundle shared by the paper's first author
  (`high_conf_homodimers_manuscript_v1.tar.gz`). The released metadata carries no quality flag
  or pLDDT, so the paper's selection cannot be re-derived; the bundle is the set. Moved from
  the extraction, with 35,254 `AF_<n>` names normalised to `AF-<n>`.
- **Heterodimers** are the recalibrated AFDB release (80,248), copied from
  `AFDB_heterodimer/cif`. Those files were gzip streams named `.cif`; here they are `.cif.gz`,
  since the CIF reader picks its decoder by extension. Bytes unchanged.
- **MSAs** are the release's own, one per monomer entity, keyed by entity number. A homodimer
  uses its own entity's MSA; a heterodimer chain uses its monomer's, linked by UniProt
  (158,267 chains) or, where the metadata has no row for the monomer, by exact sequence
  (2,229). Directory-stored files are hard links into the download mirror; tar-only ones were
  extracted. The release's `pandemic_prep/` batch is not staged — nothing here needs it.

Staged by `scripts/maintenance/stage_afdb_multimer.py` and `stage_afdb_multimer_msa.py`.

## `BioMol_clean/materials/` — derived inputs

```
materials/
├── raw/fasta/              one consolidated FASTA per DB (docs/seq-id-and-cluster-scheme.md)
│   ├── cif_pdb.fasta           PDB        extracted from pdb/cif
│   ├── ofd.fasta               OFD        short + long + rna + disordered, concatenated
│   ├── ofd_long_recovered.fasta OFD       the 14 long entries recovered 2026-09-14 (docs/long-cif-recovery-2026-09-14.md)
│   ├── teddymer.fasta          TDM        extracted from teddymer/cif
│   ├── afm.fasta               AFM        afm_homodimer.fasta + afm_heterodimer.fasta (each extracted from its cif DB)
│   ├── pdb_polypeptide_L.fasta            template search target (PDB polypeptide(L) chains)
│   └── train.fasta · valid_1.fasta        validation-stage-2 dedup inputs
└── intermediate/
    ├── msa/P/<last 3>/<seq_id>/<seq_id>.a3m      teddymer MSAs (MMseqs2, uniref30_2302)
    ├── msa_wo_lower_TDM/                         lowercase-stripped a3m, hmmbuild input
    ├── hmm_TDM/ · hmm_output_TDM/                template search (Phase 0)
    └── mmseqs_tmp/                               per-chunk scratch, emptied on success
```

`intermediate/msa` holds only teddymer's MSAs; PDB's live in
`BioMol/materials/intermediate/msa` (PDB's older intermediates -- msa, msa_wo_lower, hmm,
hmm_output -- all predate this tree and stay under `BioMol/materials/intermediate/`). Anything that scans a directory to find its queries
(the MSA and template searches do) is scoped by which tree it points at, so keep the two
apart.

## Conventions

- **Shard on the key.** Trees fan out on a slice of the entry id (the last three digits of a
  number, or the PDB id's middle two characters), so a path is computable from an id.
- **Stage, don't point at scratch.** Builds read from these trees, never from a personal
  download directory, so a scratch copy can disappear without breaking a rebuild.
- **Copy what you don't own, move what you do.** Someone else's download is copied;
  your own is moved (same filesystem, so a rename). Hard links are fine within `/data`.
- **Record the selection.** When a tree holds a subset, its README says which and why —
  the source release will not say it for you.
