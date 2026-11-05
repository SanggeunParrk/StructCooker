# Sources, sequences, ids, antigen MSA and templates

The antibody-antigen distillation set (`abag`) is built from three sources of antibody-antigen pairs, **ASD v2.0**, **AbRank** and **POISK**, with AntiBinder (light chains), SAbDab and the PDB (experimental
structures) as helpers. It is outside the `build-all` DAG, like the [interface DB](../interface-db.md): its stages live in `db_configs/abag/` and run one by one with `structcooker build abag/<name>`
([RECIPE.md](RECIPE.md)). Outputs go to their own tree, `$BIOMOL_ROOT/abag/`, so nothing of the release is rewritten.

## Sources  (`materials/raw/`, each folder has a `README.md` and a `SOURCE.tsv` with checksums)

| Source | Rows | What it gives |
|---|---:|---|
| ASD v2.0 (`asd/v2.0/`) | 1,177,819 | NaturalAntibody's antibody-antigen database: 29 datasets (buzz, alphaseq, abbd, patents, ...), a Delta table. Antigens are UniProt full chains; Kd, scores, labels, ... mixed per `affinity_type` |
| AbRank (`abrank/`) | 342,362 | a merge of nine sources (RBD-escape, CATNAP, AlphaSeq, AbCoV, SKEMPIv2, AbSci, SabDab, OVA-binders, AIntibody); a tab-separated table; the unit is in the column header |
| POISK (`poisk/zenodo_19883049/`) | 89,074 | antibody-antigen pairs an LLM pipeline extracted from patents; free-text affinity |
| AntiBinder (`antibinder/repo/`) | | ASD's covid-19, hiv, biomap and met rows as they were before ASD kept the heavy chain only: the light chains |
| SAbDab (`sabdab/`), PDB (`materials/intermediate/fasta/cif_pdb.fasta`) | | the antibody-antigen complexes with an experimental structure |

ASD is read from `asd/v2.0/` (the live Delta table staged by `scripts/staging/stage_asd.py`), never from `download/extract/` (it also holds superseded parquet files).
AbRank: the 807 rows whose heavy and light columns are swapped are put right (both columns end in the other chain's FR4); an antigen of several chains (`:`, 1,141 rows) is
not a single chain and is left out. POISK: a row POISK's own audit calls `CONTRADICTED` or `NOT_FOUND` is dropped. The binding rules per source are in
[affinity-filter.md](affinity-filter.md). Licences: ASD and AbRank (CC BY-NC 4.0 by its paper) are not checked; POISK is CC BY 4.0; AntiBinder is MIT.

## Sequences and ids  (`run.sh ids`)

| Stage | Output | What |
|---|---|---|
| `abag/antigen_fasta` | `materials/intermediate/fasta/abag_antigen.fasta` | the distinct antigen sequences of every usable row of the three sources: 21,592 |
| `abag/antibody_fasta` | `.../abag_antibody.fasta` | the distinct heavy (892,102; a VHH or scFv is a heavy) and light (72,384) chains |
| `abag/seq_id_map` | `abag/metadata/seq_id_map.tsv` | `seq_id` of every sequence above: 986,071 |

- **FASTA.** Header `ABAG<H|L|T>_<sha1[:12]> | polypeptide(L) | Auth:<H|L|A>`; the order is (role, length, sequence), because a seq_id is a counter assigned in order of
  first appearance. A usable row has a heavy chain and an antigen, every chain A-Z only; a sequence several sources hold is one record.
- **seq_id.** One id space with the release: the seed is the release's `metadata/seq_id_map.tsv` **plus** `abag/seed/seq_id_map_asd.tsv` (the ids an earlier, ASD-only run handed
  out), so a sequence the release holds keeps its release id, an ASD sequence its ASD id, and a new one continues the counter above both. None of the 880,816 ASD ids moved. The
  output holds only these sequences. When the release map is next rebuilt, `abag_antigen.fasta` and `abag_antibody.fasta` go **last** in `fasta_paths`.
- **Clusters** are made later, on the chains of the kept pairs only ([RECIPE.md](RECIPE.md) step 4); there is no clustering of every sequence.
- DB codes `ABG` (antigens) and `ABA` (antibodies), ids `cABG_<seq_id>` / `cABA_<seq_id>`.

An earlier, ASD-only run (stages `asd/...`, clusters of every ASD sequence, an ASD-only final table) was replaced by this one; its outputs are kept in `abag/legacy_asd/` and are not used.

## Antigen MSA and templates  (`run.sh antigen`)

Antigens first, antibodies later (after the pairs are filtered). The antigens are those of the kept pairs (`abag/ag_worklist`), 10-3000 aa. MSA and templates follow teddymer's path
(`db_configs/teddymer/`): MMseqs2 against UniRef30 2302 alone for the MSA, hmmsearch against the PDB polypeptide(L) chains for the templates.

| Stage | Output | What |
|---|---|---|
| `abag/ag_worklist` | `abag/worklists/seq_id_map_AG.tsv` | the antigens of the kept pairs |
| **MSA** `abag/ag_msa_search` | `abag/msa/ag/<shard>/<seq_id>.a3m` | ColabFold-style MMseqs2 search, chunks of 300, 8 nodes; an a3m that exists is skipped, a failed chunk is retried in four parts |
| `abag/ag_msa` | `abag/lmdb/msa/msa_ag.lmdb` | schema E, keyed by seq_id |
| `abag/ag_msa_wo_lower` | `abag/intermediate/msa_wo_lower_AG/` | a3m without the lowercase columns |
| `abag/ag_hmmsearch` | `abag/intermediate/hmm_output_AG/` | hmmbuild + hmmsearch against all PDB chains, no date cut (resumable) |
| **Templates, two versions** `abag/ag_template_candidates_20210930`, `_2026` | `abag/worklists/template_AG_<version>/` | up to 60 hits released up to **2021-09-30** (AF3's cutoff) and up to **2026-02-18** (the newest release of the PDB snapshot), from the same hmmsearch output |
| `abag/ag_template_20210930`, `abag/ag_template_2026` | `abag/lmdb/template/template_ag_<version>.lmdb` | schema D, keyed by seq_id, up to 20 template mols each |

Then the antigen keep list (`scripts/analysis/abag_ag_keep.py`): **an antigen is dropped when its MSA depth is 50 or less and it has no template in the 2021-09-30 version** ->
`abag/worklists/ag_keep.txt`. The final pair table only uses these antigens.

- **Polyproteins out.** An antigen over 3,000 aa is not in the table (the polyproteins and the giant single chains); peptides under 10 aa have nothing to align.
- **Searched even when the PDB has the MSA.** An antigen whose sequence the PDB holds keeps its seq_id but gets this search's MSA, so every antigen is one lineage.

## Antibody clustering on CDRs

The clusters of the final antibodies (`abag/final_antibody_seq_cluster_cdr90`, `..._cdr100`): the role (heavy / light) comes from the fasta header, ANARCI numbers each chain (Chothia, in
batches of 5,000), the CDR-H3 (heavy) and CDR-L3 (light) are clustered with CD-HIT at 0.9 or at 1.0 (identical CDRs); a chain ANARCI cannot number is clustered with mmseqs2 at 0.40.
At 1.0 `run_CDHIT` also sets cd-hit's length and coverage cutoffs to 1.0: its identity is over the shorter sequence, so with the 0.8 cutoffs a CDR that is a piece of a longer one
would join it. Ids are `cABA_<rep seq_id>`. An antigen cluster at 100% is the seq_id itself, so "antigen 100% x CDR 100%" needs no antigen clustering.
