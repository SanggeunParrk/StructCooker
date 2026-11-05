# Recipe: the antibody-antigen distillation pairs

How the pair tables `abag/pairs/pairs_basic.tsv` and `abag/final/pairs_final.tsv` are made from raw data, stage by stage, with the command, the inputs, the output and what to check.
Every stage is a DataCooker recipe: a `db_configs/abag/<name>.yaml` plus a recipe in `src/structcooker/workflows/`, built with
`pixi run structcooker build abag/<name>` (an sbatch job; nothing heavy runs on the login node). `run.sh` in this folder chains them.
The reasons behind the choices are in [pair-filters.md](pair-filters.md), [affinity-filter.md](affinity-filter.md),
[light-chain-pairing.md](light-chain-pairing.md) and [sources.md](sources.md); [overview.html](overview.html) shows the whole thing.

Three sources: **ASD v2.0, AbRank, POISK**. All paths are under `$BIOMOL_ROOT` (default `/data/shared/cssb_data/BioMol_clean`).

## 0. Raw data  (`materials/raw/`)

| Source | How | Result |
|---|---|---|
| ASD v2.0 | the Drive folder `1iSF3BZdDu2jjxIz3sSz5wjJMbPsAkPDM`, **downloaded by hand in a browser** (gdown fails on the quota); put the four zips in `asd/download/`, then `sbatch scripts/staging/unpack_asd.sbatch` and `sbatch scripts/staging/stage_asd.sbatch` | `asd/v2.0/` (the live Delta table, 20 parquet files, 1,177,819 rows), `asd/SOURCE.tsv` |
| AbRank | `sbatch scripts/staging/copy_abrank.sbatch` (copy from cssb-master, sha256 checked; the Kaggle original needs a login) | `abrank/AbRank_dataset_updated_v2.csv` (a TSV, 342,362 rows) |
| POISK | `sbatch scripts/staging/download_poisk.sbatch` (Zenodo 19883049, md5 checked) | `poisk/zenodo_19883049/poisk_dataset.csv` (89,074 rows) |
| AntiBinder | `sbatch scripts/staging/download_antibinder.sbatch` (git clone, MIT) | `antibinder/repo/` |
| SAbDab summary | `sabdab/sabdab_summary_all.tsv` (already there); the PDB sequences are `materials/intermediate/fasta/cif_pdb.fasta` | |

Each folder has a `README.md` (provenance, licence) and a `SOURCE.tsv` (checksums). ASD and AbRank licences are not checked (AbRank is CC BY-NC 4.0 by its paper).
The seed of the id space `abag/seed/seq_id_map_asd.tsv` is the ASD ids of the earlier ASD-only run (a copy of its `seq_id_map.tsv`): it must exist before step 1.

## 1. Sequences and ids  (`run.sh ids`, about 5 minutes)

| Stage | Output | Notes |
|---|---|---|
| `abag/antigen_fasta` | `materials/intermediate/fasta/abag_antigen.fasta` | the distinct antigens of every usable row of the three sources (21,592) |
| `abag/antibody_fasta` | `.../abag_antibody.fasta` | distinct heavy (892,102) and light (72,384) chains; a VHH or scFv is a heavy |
| `abag/seq_id_map` | `abag/metadata/seq_id_map.tsv` (986,071 ids) | one id space with the release; the seed is the release map plus the ASD ids, so no id moves (checked: none of 880,816 ASD ids changed) |

## 2. The pairs  (`run.sh pairs`, about 5 minutes)

| Stage | Output | What it does |
|---|---|---|
| `abag/pairs_binding` | `abag/pairs/pairs_binding.tsv` (677,817 pairs) | S2 (heavy chain and antigen present, A-Z), antigen 10-3000 aa, S3 (distinct heavy/light/antigen sequence triples across the three sources), the binding filter Kd <= 1 uM per source ([affinity-filter.md](affinity-filter.md)). **No S1**; no MSA filter yet |
| `abag/pairs_chains` | `abag/pairs/pairs_chains.tsv` | chain status (paired / scfv / vhh / filled from AntiBinder / ambiguous / unresolved) and the recovered light sequence |
| `abag/pairs_basic` | `abag/pairs/pairs_basic.tsv` (676,517 rows) | duplicates merged, polyprotein, PDB structure (heavy chain, light chain and antigen all identical to one SAbDab entry), and the `keep` flag |

`keep == 1`: **607,342 pairs** (ASD 484,207, AbRank 122,984, POISK 10,593; a pair two sources report counts in both), 14,751 antigens. Columns of `pairs_basic.tsv`:
`antigen_seq_id, heavy_seq_id, light_seq_id, light_sequence, chain_status, light_source, datasets (source:dataset), polyprotein, pdb_structure, keep`.
A different choice of filter is a different selection on these columns, not a rebuild.

## 3. Antigen MSA and templates, two versions  (`run.sh antigen`, many hours)

Antigens first; the antibodies' MSA and templates come after the filtering. The antigens are those of the kept pairs; the 5,998 that were searched in the earlier ASD-only run are skipped
by the MSA search (an a3m that exists is kept), the other 8,753 (median 223 aa) are searched.

| Stage | Output |
|---|---|
| `abag/ag_worklist` | `abag/worklists/seq_id_map_AG.tsv`: the antigens of the kept pairs |
| `abag/ag_msa_search` | `abag/msa/ag/<shard>/<seq_id>.a3m`: MMseqs2 against UniRef30 2302 only, ColabFold style, chunks of 300 on 8 nodes |
| `abag/ag_msa` | `abag/lmdb/msa/msa_ag.lmdb` |
| `abag/ag_msa_wo_lower` -> `abag/ag_hmmsearch` | `abag/intermediate/hmm_output_AG/`: hmmsearch against all PDB chains, no date cut |
| `abag/ag_template_candidates_20210930` and `..._2026` | per seq_id candidates released up to **2021-09-30** and up to **2026-02-18**, 60 candidates, from the same hmmsearch output |
| `abag/ag_template_20210930` and `abag/ag_template_2026` | `abag/lmdb/template/template_ag_20210930.lmdb`, `template_ag_2026.lmdb`: schema D, up to 20 templates per seq_id |

Then the keep list of antigens (`scripts/analysis/asd_msa_depth.py`, `asd_template_counts.py`, `abag_ag_keep.py`): **an antigen is dropped when its MSA depth is 50 or less and it has no
template in the 2021-09-30 version** -> `abag/worklists/ag_keep.txt`. The 2026 templates are for use; the filter is on the 2021 ones so that it does not depend on structures newer than the AF3 cutoff.

## 4. The final chains: fasta, seq_ids, clusters, the final pair table  (`run.sh final`, about 40 minutes)

Steps 1-2 numbered every sequence; this step does it again for the chains of the **kept pairs whose antigen is on the keep list** (`pairs_basic.tsv`, `keep == 1`, `ag_keep.txt`), recovered
light chains included, and ends in one table with a seq_id and a cluster for every chain. Outputs are under `abag/final/`.

| Stage | Output | What |
|---|---|---|
| `abag/final_antigen_fasta` | `abag/final/fasta/antigen.fasta` | the antigens of these pairs |
| `abag/final_antibody_fasta` | `abag/final/fasta/antibody.fasta` | their heavy and light chains, the recovered light chains included |
| `abag/final_seq_id_map` | `abag/final/metadata/seq_id_map.tsv` | seeded with `abag/metadata/seq_id_map.tsv`: a chain keeps its id, a recovered light chain continues the counter |
| `abag/final_antigen_seq_cluster40` | `abag/final/metadata/seq_cluster40_ABG.tsv` | mmseqs2, 40% identity, coverage 0.8, on the final antigens |
| `abag/final_antibody_seq_cluster_cdr90` | `.../seq_cluster_cdr90_ABA.tsv` | CDR-H3 / CDR-L3 (ANARCI + CD-HIT) at 0.9 |
| `abag/final_antibody_seq_cluster_cdr100` | `.../seq_cluster_cdr100_ABA.tsv` | the same at 1.0 (identical CDRs) |
| `abag/final_pairs` | `abag/final/pairs_final.tsv` | the pairs with seq_ids and clusters (`antigen_cluster40`, `heavy_cluster_cdr100`, `light_cluster_cdr100`, `heavy_cluster_cdr90`, `light_cluster_cdr90`) |

The two fastas first, then `final_seq_id_map`, then the three clusterings (in parallel), then `final_pairs`. Cluster ids are `cABG_<rep seq_id>` / `cABA_<rep seq_id>` **of this subset**.

## 5. Checks and analyses  (`scripts/analysis/`)

| Script | Answers |
|---|---|
| `asd_final_stats.py` | counts of the final pair table: chains, clusters, antigen x heavy cluster combinations and their sizes |
| `asd_pairing_gold.py` | heavy-light pairing against SAbDab: error rate of ASD patents, POISK, AntiBinder |
| `asd_structure_rules.py` | the structure rule as exact identity, as containment, and ASD's `target_pdb`, side by side |
| `asd_dataset_funnel.py`, `asd_antigen_funnel.py` | where each ASD dataset / antigen is lost on the way |
| `asd_*`, one-off | questions answered during the work (ASD only) |

Tests: `pixi run pytest tests/instructions` (the readers `abag.py`, `abag_binding.py`, `abag_pairs.py`, `abag_chains.py`, `abag_basic.py`, `abag_final.py`, `asd_binding.py`).

## Reproducing and changing it

- `run.sh <group>` rebuilds a group from the outputs of the groups before it.
- **Change a binding rule**: `src/structcooker/instructions/readers/asd_binding.py` (ASD) or `abag_binding.py` (AbRank, POISK), then rebuild from `abag/pairs_binding`.
- **Keep the structure-carrying pairs**: use `keep` without the `pdb_structure` condition; **a looser structure test**: `structure_match: contained` in `db_configs/abag/pairs_basic.yaml`.
- **Fill light chains from POISK too**: a new source in `abag_chains.py` (not done: about 14% of its unique candidates disagree with SAbDab).
- **A different template cutoff**: copy `db_configs/abag/ag_template_candidates_2026.yaml` and `ag_template_2026.yaml` with another `date_cutoff`.

## Not done

A per-cluster cap of the final table (the libraries: abbd, alphaseq, buzz), the antibody MSA and templates (and the MSA of the recovered light chains), the licences, and an ipTM comparison of the
filled pairs.
