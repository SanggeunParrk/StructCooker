# End-to-end build

From downloaded sources to every BioMol_clean DB: PDB, OpenFold distillation (OFD),
teddymer (TDM) and AFDB multimer (AFM). Downloading is described, not automated here;
everything after it is `structcooker`.

```
download (manual, §1)  ->  stage (scripts/staging, §2)  ->  reference inputs in place (§3)
    ->  structcooker build-all --manifest db/MANIFEST_all.yaml  (§4)
```

`DATA_ROOT` defaults to `/data/shared/cssb_data`, `OUTPUT_ROOT` to `DATA_ROOT/BioMol_clean`.
The trees and the raw/intermediate rule are in [raw-materials.md](raw-materials.md):
**`BioMol/materials/raw/` holds exactly what was downloaded; anything we derive goes to an
`intermediate/` tree.**

## 1. Download

| Source | Get it from | Put it at (`DATA_ROOT/BioMol/materials/raw/…`) |
|---|---|---|
| PDB mmCIF | `structcooker download mmcif` (wwPDB rsync) | `cif/<mid 2>/<id>.cif.gz` |
| CCD | `structcooker download ccd` | `ccd/components_<date>.cif.gz` (the split goes to `intermediate/ccd/components/`) |
| SabDab summary | `structcooker download sabdab` | `DATA_ROOT/external/SabDab/` |
| OpenFold distillation | `structcooker download openfold`, or `s5cmd --no-sign-request` from `s3://openfold3-data` | `openfold_distillation/` (keep the S3 layout) |
| ↳ OFD short/long relaxed models | per id: `s5cmd cp 's3://openfold3-data/monomer_distillation_sets_v2/{short,long}_monomers/raw/<id>/best_structure_relaxed.pdb*' <dst>/<id>/` — only this file is needed (≈5.8 TB for long; the whole long `raw/` is ≈140 TB). 80 long ids have no model in the release; they are left out. | `openfold_distillation/monomer_distillation_sets_v2/{short,long}_monomers/raw/<id>/` |
| ↳ RNA alignment arrays | `s3://openfold3-data/pdb_training_set/alignment_arrays/<pdb>_<chain>.npz` | `rna_alignment_arrays/` |
| Teddymer | `https://teddymer.steineggerlab.workers.dev/foldseek/teddymer_dimerpdbs.tar` + the metadata tarball (`cluster.tsv`, `nonsingletonrep_metadata.tsv`) | a scratch dir; §2 stages the selection |
| AFDB multimer | `https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/` — `homodimer_metadata.csv`, `heterodimer_metadata.csv`, the heterodimer models (`passes_quality_threshold`), `msas/` (batches 251124_1.2M, 251124_8.8M, 251209_13.4M). Homodimers: the manuscript bundle `high_conf_homodimers_manuscript_v1.tar.gz` from the paper's authors. | a scratch dir; §2 stages them |

OpenFold publishes no pLDDT for the RNA and disordered sets (checked 2026-10-02: npz, cif,
training caches and portal indexes), so those DBs carry b_factor 0.

## 2. Stage (scripts/staging)

Run as batch jobs. Each writes only to `raw/` (selection, layout, names) or to `intermediate/`
(anything it computes).

| Step | Command | Result |
|---|---|---|
| Teddymer selection | `python scripts/staging/stage_teddymer.py --teddymer-root <dl> --out $DATA_ROOT/BioMol/materials/raw/teddymer` | 510,454 non-singleton representatives passing the source paper's interface filter, `pdb/<last 3>/…`, `SOURCE.tsv`, `README.md` |
| AFM structures | `python scripts/staging/stage_afdb_multimer.py --out $DATA_ROOT/BioMol/materials/raw/afdb_multimer` | homodimer `pdb/`, heterodimer `cif/` (`.cif.gz`), `AF_`→`AF-` names, `SOURCE.tsv` |
| AFM MSAs | `stage_afdb_multimer_msa.py plan`, then `stage --task i --ntasks N` (array), then `finalize` (`AFM_MSA_DOWNLOAD`, `AFM_STAGE_WORK`, `AFM_HET_METADATA` point at the download) | release MSAs in `raw/afdb_multimer/msa/`; `intermediate/afdb_multimer/heterodimer/chain_msa.tsv` (content-verified chain → MSA) |
| AFM missing MSA | `sbatch scripts/staging/generate_afm_missing_msa.sbatch AF-<n>` for each entity the release lacks (`plan` lists them; 1 in 2026-09) | `intermediate/afdb_multimer/msa_generated/AF-<n>-msa_v1.a3m[.zst]`; `finalize` above picks it up |
| AFM mmCIF repairs | `python scripts/staging/repair_afdb_entity_poly.py` | 6 heterodimers whose `_entity_poly` contradicts `_entity_poly_seq`: repaired copies + `overrides.json` under `intermediate/afdb_multimer/heterodimer/repaired/` (raw untouched) |

## 3. Reference inputs

Supplied, not derivable from the downloads; each fixes part of the result to an earlier one.
Without them a build is self-consistent but not identical to this release.

| Input | Location | Fixes |
|---|---|---|
| seq_id seed | `DATA_ROOT/reference/seq_id_map.tsv` | production seq_ids; new sequences are appended in `db/metadata/seq_id_map.yaml` order (PDB, OFD, TDM, AFM, recovered long) |
| OFD reference fastas | `DATA_ROOT/BioMol/materials/intermediate/openfold_distillation/fasta/{short,long,rna,disordered}.fasta` | the OFD sequence set seq_id_map and OFD clustering run over |
| PDB MSAs | `DATA_ROOT/BioMol/materials/intermediate/msa/` | PDB a3m DBs and the PDB template search |
| SignalP results | `OUTPUT_ROOT/materials/intermediate/signalp/` | sequence trimming |
| Historical CCD | `OUTPUT_ROOT/metadata/reference_ccd.lmdb` | 24 component definitions older PDB entries need |
| Chain model choice | `OUTPUT_ROOT/metadata/reference_chain_selection.lmdb` | production's per-chain model choice on occupancy ties |
| Disordered set work lists | `DATA_ROOT/BioMol/materials/intermediate/openfold_distillation/disordered_set/` | which disordered MSAs / templates are built |

## 4. Build

```bash
structcooker inspect --strict --manifest db/MANIFEST_all.yaml     # every input present, DAG closed
structcooker build-all --manifest db/MANIFEST_all.yaml --workdir <scratch>
```

`db/MANIFEST_all.yaml` is the union of the per-set manifests (`MANIFEST.yaml`,
`MANIFEST_{pdb,distillation,teddymer,afdb_multimer,msa,cifcore*}.yaml`); regenerate it after
editing one (`python -c 'from structcooker.manifests import write; write()'`, enforced by a
test). Each set's manifest still builds that set alone.

Sizes and times on this cluster (2026-09/10): OFD long cif ≈ 8 × 3 h; teddymer MSA (MMseqs2,
8 nodes) ≈ 2 days; teddymer hmmsearch (8 nodes) ≈ 31 h; AFM hmmsearch (8 nodes) ≈ 6 days.
Large `parallel` searches declare `node_count: 8` and are submitted as an 8-task array; each
task refuses to start on a node whose `/dev/shm` cannot hold process semaphores (joblib would
otherwise run serially there).
