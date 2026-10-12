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
| OFD long npz recovery | `python scripts/staging/recover_ofd_long_npz.py` | 14 long `structure.npz` that end in stale bytes (2026-09-14): their valid zip prefix, hash-checked, + `overrides.json` under `intermediate/openfold_distillation/long_recovered/` (raw untouched); read by the long template build |
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

A stage may fail only on inputs listed in `db/release-exceptions.json`, each with its reason
and SHA256: 5 PDB entries the parser cannot map, and 840 OpenFold `alignment.npz` files that
are empty or corrupt in the release (OFD long 667, short 173; the MSA DBs lack those
sequences). Any other failure stops the stage.

`db/MANIFEST_all.yaml` is the union of the per-set manifests (`MANIFEST.yaml`,
`MANIFEST_{pdb,distillation,teddymer,afdb_multimer,msa,cifcore*}.yaml`); regenerate it after
editing one (`python -c 'from structcooker.manifests import write; write()'`, enforced by a
test). Each set's manifest still builds that set alone.

Sizes and times on this cluster (2026-09/10): OFD long cif ≈ 8 × 3 h; teddymer MSA (MMseqs2,
8 nodes) ≈ 2 days; teddymer hmmsearch (8 nodes) ≈ 31 h; AFM hmmsearch (8 nodes) ≈ 6 days.
Large `parallel` searches declare `node_count: 8` and are submitted as an 8-task array; each
task refuses to start on a node whose `/dev/shm` cannot hold process semaphores (joblib would
otherwise run serially there); set `DATACOOKER_EXCLUDE_NODES=node04,node07` (any such
nodes) before `build-all` so the scheduler does not place tasks on them.

## 5. Sampled end-to-end check

A full build is days and terabytes. To check that the release DAG still builds end to
end after a change, run every stage on a sample of a few hundred inputs instead
(`structcooker.sample_e2e`). All four commands are login-node safe: `prepare` and `verify`
submit SLURM jobs, and `run` submits stage jobs the way `build-all` does.

```bash
R=/data/shared/cssb_data/e2e_sample/$(date +%Y%m%d)      # any new, empty directory
structcooker e2e-sample prepare $R      # job: pick the sample, build $R/data and $R/out
structcooker e2e-sample run $R          # build-all on the sample, run state in $R/run
structcooker e2e-sample status $R       # stage phases; rerun `run` to resume after a fix
structcooker e2e-sample verify $R       # job: compare with the release -> $R/verify.json
```

**Sample.** By default 30 PDB entries and 20 entries each from OFD short, long, rna and
disordered, teddymer, and AFM homo- and heterodimers (`--count pdb=50` etc.). The
selection is a seeded hash order, so the same production data gives the same sample.
On top of that random selection, the sample also includes these cases on purpose:

- the accepted PDB failures (`db/pdb-exceptions.json`)
- two manually fixed mmCIFs
- two entries from each published PDB-derived DB (train, valid1/2, mpnn, chain), so every filter has something to keep
- the six repaired AFM heterodimers
- the homodimer whose MSA we generated
- two recovered OFD long entries

Their files are linked one by one into `$R/data` (DATA_ROOT) in the production layout.
The reference fastas, disordered work lists, `chain_msa.tsv` and `overrides.json` are
rewritten for the sample. The following are linked whole:

- the search DBs
- the CCD
- the seq_id seed
- the reference CCD and chain-selection DBs

**Run.** `DATA_ROOT`/`OUTPUT_ROOT` point at the sample, every other root override is
cleared, and `run` refuses to start if any config still resolves outside `$R`. Parallel
searches run on one node (`STRUCTCOOKER_MAX_NODES=1`). Nothing in production is written.

**Verify.** For every stage `verify.json` records whether it completed and, per output,
how it compares with the release at the same key:

- **Must be identical (fails otherwise).** Records that depend only on their own inputs:
  structures, MSAs and the CCD (schemas A, E and G).
- **Reported only.** Records that depend on the rest of the build. Examples: attached
  metadata (clusters are computed over the sample), templates (searched against the
  sample's PDB chains), and train/valid sets. Their differing fields are counted, so an
  unexpected one stands out.
- **Projections.** These are compared line by line.

A run is pinned to the code it started with (the release definition includes a snapshot of
`src/`). After a fix, `run` on the same root refuses ("Release definition changed"):
`prepare` a new root (about 5 minutes) and run that. A stage shown as `stalled` lost its
work job; its attempt logs are under `$R/run/nodes/<stage>/attempts/`.
