# Raw materials — what the builds read, and where it came from

Two trees hold build inputs, and they are not the same kind of thing:

| tree | holds | written by |
|---|---|---|
| `DATA_ROOT/BioMol/materials/raw/` | **exactly what was downloaded** (a selected subset may be copied; contents are never edited) | staging scripts, once |
| `DATA_ROOT/BioMol/materials/intermediate/` | anything **we made** from those sources: extracted fastas, repaired files, maps, generated MSAs, legacy search outputs | staging / maintenance scripts |
| `OUTPUT_ROOT/materials/intermediate/` (`BioMol_clean`) | what the pipeline makes: per-DB FASTAs, MSAs, search outputs | the pipeline |

Every staged source directory carries a `README.md` (what it is, why this subset) and a
`SOURCE.tsv` (where each file came from, its size). The layout below is the contract the
`db/*.yaml` configs rely on; a README in the directory itself is the detailed record.

## `BioMol/materials/raw/` — source data

```
raw/
├── ccd/                    components_<date>.cif.gz              wwPDB CCD snapshot (as downloaded)
├── cif/                    <mid 2 chars>/<pdb id>.cif.gz         RCSB mmCIF mirror, 1,101 shards
├── rna_alignment_arrays/   <pdb>_<chain>.npz · DOWNLOAD.log      from s3://openfold3-data/pdb_training_set
│
├── openfold_distillation/                                        OFD — OpenFold3-preview2 distillation
│   ├── monomer_distillation_sets_v2/   README.md (OpenFold's) · shared/reference_mols/
│   │   ├── short_monomers/   preprocessed/ 430,420 · raw/ 430,420 · cache.json · cache_lmdb/
│   │   └── long_monomers/    preprocessed/ 16,099,486 · raw/ 8,083 · cache.json · cache_lmdb/
│   ├── rna_distillation_set/            rna_monomer_preprocessed_cache/ 126,780 (no raw)
│   └── disordered_set/                  structure_files/ · templates/ 28,569 · alignment_arrays/
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
    │   └── cif/<last 3>/AF-<n>-model_v1.cif.gz
    └── msa/                1,799,836  monomer MSAs the release ships, for the entities used
    │   ├── README.md
    │   ├── SOURCE.tsv      entity · release batch dir / tar member / generated · used by
    │   └── <last 3>/AF-<n>-msa_v1.a3m.zst
```

### teddymer

Selected, not the full release: 510,454 non-singleton cluster representatives that pass the
source paper's interface filter (length > 10, pAE < 10, pLDDT > 70) out of 10,089,503 dimers.
Each file is two TED domains of ONE AFDB v4 model — a domain-domain interface, not two
proteins. A domain can be sequence-discontinuous, so a chain's residue numbers may jump; the
jump is a domain boundary, not missing structure. B-factor holds pLDDT.
Staged by `scripts/maintenance/stage_teddymer.py` (copy).

### openfold_distillation

The OpenFold3-preview2 distillation release, downloaded from its S3 bucket (`bin/s5cmd`).
Builds read `preprocessed/` (`structure.npz` / `alignment.npz` / `template.npz` per entry).
`raw/` is the pipeline's own output (`best_structure_relaxed.pdb`, with pLDDT in the B column
and hydrogens; MSAs; `hmm_output.sto`), complete for short but only 8,083 entries for long.
`structure.npz` carries no B-factor or pLDDT, so the OFD cif DBs hold b_factor = 0.
Moved here 2026-09-30 from `/data/shared/cssb_data/openfold_distillation`, which is now a
link to this directory so older paths keep resolving.

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

## `BioMol/materials/intermediate/` — made from the sources

```
intermediate/
├── fasta/                    cif.fasta · cif_pdb.fasta · pdb_polypeptide_L.fasta (extracted from PDB)
├── ccd/components/           <COMP_ID>.cif × 50,782, split from raw/ccd/components_20260803.cif.gz
├── afdb_multimer/            README.md
│   ├── heterodimer/chain_msa.tsv                 chain -> monomer MSA entity
│   ├── heterodimer/repaired/                     6 repaired mmCIFs + ENTITY_POLY_REPAIRS.tsv + overrides.json
│   └── msa_generated/                            the 1 MSA the release does not ship
├── openfold_distillation/    README.md · fasta/ · disordered_set/ lists · lmdb/ · bin/ · _archive/
└── msa/ · msa_wo_lower/ · hmm/ · hmm_output/ · hhm/ · hhr/ · signalp/ · …   PDB's legacy search outputs
```

A build that needs both a release file and our change reads the raw file and applies the
change from here -- AFM heterodimers through `AFM_CIF_OVERRIDES`, the AFM MSA DB through
`metadata/afm_msa_filelist` (raw MSAs + the generated one) -- so raw is never edited.

## `BioMol_clean/materials/` — what the pipeline makes

```
materials/intermediate/
├── fasta/                  one consolidated FASTA per DB (docs/seq-id-and-cluster-scheme.md)
│   ├── cif_pdb.fasta           PDB        extracted from pdb/cif
│   ├── ofd.fasta               OFD        short + long + rna + disordered, concatenated
│   ├── ofd_long_recovered.fasta OFD       the 14 long entries recovered 2026-09-14
│   ├── teddymer.fasta          TDM        extracted from teddymer/cif
│   ├── afm.fasta               AFM        afm_homodimer.fasta + afm_heterodimer.fasta
│   ├── pdb_polypeptide_L.fasta            template search target (PDB polypeptide(L) chains)
│   └── train.fasta · valid_1.fasta        validation-stage-2 dedup inputs
├── msa/P/<last 3>/<seq_id>/<seq_id>.a3m   teddymer MSAs (MMseqs2, uniref30_2302)
├── msa_wo_lower_{TDM,AFM}/                lowercase-stripped a3m, hmmbuild input
├── hmm_{TDM,AFM}/ · hmm_output_{TDM,AFM}/ template search (Phase 0)
├── hmm_reduced_TDM/                       Phase 1/2 reduced hits
└── mmseqs_tmp/                            per-chunk scratch, emptied on success
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
- **raw is what was downloaded, nothing else.** Extracting a tar, sharding, renaming and
  choosing a subset are fine; editing contents, adding files we generated, or keeping our
  maps and lists there is not -- those go to `intermediate/`, and a build that needs a change
  applies it at read time (overrides, file lists). Exception still to fix: `structcooker
  fix-cif` writes the manual PDB fixes into the mmCIF tree itself.
