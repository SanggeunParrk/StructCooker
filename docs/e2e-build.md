# End-to-end build

From downloaded sources to every BioMol_clean DB: PDB, OpenFold distillation (OFD),
teddymer (TDM) and AFDB multimer (AFM). Part of the downloading is `structcooker download`,
the rest is described here; everything after it is `structcooker`. All of it runs as SLURM
jobs; the commands below only submit them, so they are safe on a login node.

```
download (§1)  ->  stage (scripts/staging, §2)  ->  reference inputs in place (§3)
    ->  structcooker build-all --manifest db/MANIFEST_all.yaml --workdir RUN  (§4)
    ->  release-check, release-docs (§5)
```

§6 checks the DAG on a sample instead of the full data; §7 fixes some records of a published
DB.

Everything lives under one root, `BIOMOL_ROOT` (default `/data/shared/cssb_data/BioMol_clean`):
inputs under `materials/`, DBs under `lmdb/`, released tables under `metadata/`, the release's
own docs under `docs/`. Only the shared UniRef30/MMseqs2 search databases sit outside it,
under `SEQ_DB_ROOT` (default `/data/shared/cssb_data/seq_db`); the release DAG does not read
them (only the hand-run `teddymer/msa_search` and the AFM missing-MSA job do). Cluster
settings come from the environment: `source scripts/cluster/cssb.env` (see
[development](development.md#cluster-settings)). The tree and the raw/reference/intermediate
rule are in [raw-materials.md](raw-materials.md): **`materials/raw/` holds exactly what was
downloaded; `materials/reference/` holds supplied inputs that cannot be derived from the
downloads; anything derived goes to `materials/intermediate/`.**

## 1. Download

| Source | Get it from | Put it at (`BIOMOL_ROOT/materials/raw/…`) |
|---|---|---|
| PDB mmCIF | `structcooker download mmcif` (wwPDB rsync) | `cif/<mid 2>/<id>.cif.gz` |
| CCD | `structcooker download ccd` | `ccd/components_<date>.cif.gz` (the split goes to `materials/intermediate/ccd/components/`) |
| SabDab summary | `structcooker download sabdab` | `sabdab/sabdab_summary_all.tsv` |
| OpenFold distillation | `structcooker download openfold`, or `s5cmd --no-sign-request` from `s3://openfold3-data` | `openfold_distillation/` (keep the S3 layout) |
| ↳ OFD short/long relaxed models | per id: `s5cmd cp 's3://openfold3-data/monomer_distillation_sets_v2/{short,long}_monomers/raw/<id>/best_structure_relaxed.pdb*' <dst>/<id>/` — only this file is needed (≈5.8 TB for long; the whole long `raw/` is ≈140 TB). 80 long ids have no model in the release; they are left out. | `openfold_distillation/monomer_distillation_sets_v2/{short,long}_monomers/raw/<id>/` |
| ↳ RNA alignment arrays | `s3://openfold3-data/pdb_training_set/alignment_arrays/<pdb>_<chain>.npz` | `rna_alignment_arrays/` |
| Teddymer | `https://teddymer.steineggerlab.workers.dev/foldseek/teddymer_dimerpdbs.tar` + the metadata tarball (`cluster.tsv`, `nonsingletonrep_metadata.tsv`) | a scratch dir; §2 stages the selection |
| AFDB multimer | `https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/` — `homodimer_metadata.csv`, `heterodimer_metadata.csv`, the heterodimer models (`passes_quality_threshold`), `msas/` (batches 251124_1.2M, 251124_8.8M, 251209_13.4M). Homodimers: the manuscript bundle `high_conf_homodimers_manuscript_v1.tar.gz` from the paper's authors. | a scratch dir; §2 stages them |

OpenFold publishes no pLDDT for the RNA and disordered sets (checked 2026-10-02: npz, cif,
training caches and portal indexes), so those DBs carry b_factor 0.

## 2. Stage (scripts/staging)

Run as batch jobs. The scripts default to `BIOMOL_ROOT` (`generate_afm_missing_msa.sbatch` also
reads `SEQ_DB_ROOT`). Each writes only to `materials/raw/` (selection, layout, names) or to
`materials/intermediate/` (anything it computes).

| Step | Command | Result |
|---|---|---|
| Teddymer selection | `python scripts/staging/stage_teddymer.py --teddymer-root <dl> --out $BIOMOL_ROOT/materials/raw/teddymer` | 510,454 non-singleton representatives passing the source paper's interface filter, `pdb/<last 3>/…`, `SOURCE.tsv`, `README.md` |
| AFM structures | `python scripts/staging/stage_afdb_multimer.py --out $BIOMOL_ROOT/materials/raw/afdb_multimer` | homodimer `pdb/`, heterodimer `cif/` (`.cif.gz`), `AF_`→`AF-` names, `SOURCE.tsv` |
| AFM MSAs | `stage_afdb_multimer_msa.py plan`, then `stage --task i --ntasks N` (array), then `finalize` (`AFM_MSA_DOWNLOAD`, `AFM_STAGE_WORK`, `AFM_HET_METADATA` point at the download) | release MSAs in `raw/afdb_multimer/msa/`; `intermediate/afdb_multimer/heterodimer/chain_msa.tsv` (content-verified chain → MSA) |
| AFM missing MSA | `sbatch scripts/staging/generate_afm_missing_msa.sbatch AF-<n>` for each entity the release lacks (`plan` lists them; 1 in 2026-09) | `intermediate/afdb_multimer/msa_generated/AF-<n>-msa_v1.a3m[.zst]`; `finalize` above picks it up |
| OFD long npz recovery | `python scripts/staging/recover_ofd_long_npz.py [--root $BIOMOL_ROOT]` | 14 long `structure.npz` that end in stale bytes (2026-09-14): their valid zip prefix, hash-checked, + `overrides.json` under `intermediate/openfold_distillation/long_recovered/` (raw untouched); read by the long template build |
| AFM mmCIF repairs | `python scripts/staging/repair_afdb_entity_poly.py` | 6 heterodimers whose `_entity_poly` contradicts `_entity_poly_seq`: repaired copies + `overrides.json` under `intermediate/afdb_multimer/heterodimer/repaired/` (raw untouched) |

## 3. Reference inputs

Supplied, not derivable from the downloads; each fixes part of the result to an earlier one.
They all live under `BIOMOL_ROOT/materials/reference/`. Without them a build is
self-consistent but not identical to this release.

| Input | Location | Fixes |
|---|---|---|
| seq_id seed | `seq_id_map_seed.tsv` | this release's seq_id_map (all DBs); a rebuild keeps every id, and new sequences are appended in `db/metadata/seq_id_map.yaml` order |
| OFD reference fastas | `openfold_distillation/fasta/{short,long,rna,disordered}.fasta` | the OFD sequence set seq_id_map and OFD clustering run over |
| PDB MSAs | `msa/` | PDB a3m DBs and the PDB template search |
| Teddymer MSAs | `msa_teddymer/` | teddymer MSA DBs and template search (MMseqs2, searched once; `teddymer/msa_search` is run by hand for new sequences only) |
| SignalP results | `signalp/` | sequence trimming |
| Historical CCD | `reference_ccd.lmdb` | 24 component definitions older PDB entries need |
| Chain model choice | `reference_chain_selection.lmdb` | production's per-chain model choice on occupancy ties |
| Disordered set work lists | `openfold_distillation/disordered_set/` | which disordered MSAs / templates are built (paths relative to `materials/raw/openfold_distillation`) |

## 4. Build

```bash
source scripts/cluster/cssb.env                                    # BIOMOL_ROOT, SLURM settings
structcooker inspect --strict --manifest db/MANIFEST_all.yaml      # every input present, DAG closed
structcooker build-all --manifest db/MANIFEST_all.yaml --workdir <run dir>
```

`--manifest` defaults to `db/MANIFEST_all.yaml`, the whole release; a per-set manifest builds
that set alone. `--dry-run` validates the DAG and submits nothing. The command returns once the ready
stages are submitted; there is no resident controller. `<run dir>/summary.json` shows each
stage's phase, and `complete: true` there is what marks the release built.

**What a stage does.** On its first invocation `build-all` freezes the release definition
into the run dir: a snapshot of `src/`, `libs/datacooker/src`, `pyproject.toml` and
`pixi.lock`, every stage config with its paths resolved, and `db/release-exceptions.json`.
Every job runs that snapshot. Each stage then runs as finite SLURM jobs:

1. a **plan** job fingerprints the stage's inputs (raw payloads, file lists, reader maps),
   claims its outputs and submits the work: the LMDB pipeline (shards → merge → index) for
   build / rebuild, one job or an N-task array for the workflow ops;
2. a **verify** job (a full node) decodes every output record and checks it against the
   schema, the shard accounting and the input → output key coverage;
3. only a verified stage releases its dependents, which the verify job submits.

Outputs of at most 5,000,000 records are merged into one LMDB; larger ones are published as a
shard collection (`shards.json`), see [development](development.md#lmdb-storage).

**Resume.** Rerun the same command with the same roots, code and run dir. Jobs still active
are kept; a stage that ended is replayed from its recorded attempt, rerunning only the
worker shards that did not finish, then publication and indexing; a completed stage has its
input and output fingerprints rechecked on a compute node before it is reused. If the code,
a config or the exceptions changed, the command refuses ("Release definition changed"): use a
new run dir and a new `BIOMOL_ROOT`. An output that exists but was not produced by this run
is refused too ("Untracked output"), so a release is always built into an empty tree. A
submission whose outcome cannot be established (no job id recorded) fails closed and is not
resubmitted; inspect the stage's `<run dir>/nodes/<stage>/` and its attempt logs. Do not
delete the run dir or shard directories behind a published collection.

**Failures.** A stage may fail only on inputs listed in `db/release-exceptions.json`, each by
exact key with its reason and SHA256 (verification rechecks the hash): the 5 PDB entries the
parser cannot map, and 840 OpenFold `alignment.npz` files that are empty or corrupt in the
release (OFD long 667, short 173; the MSA DBs lack those sequences). Any other failed record
stops the stage, and a rebuild may not fail any record. Two inputs with the same key fail the
stage, except where a config sets `duplicate_input_policy: first_path` (the OFD MSA and
short/long cif stages, the AFM MSA): the lexicographically first path is kept and the choice
recorded in the attempt's `input-selection.json`.

`db/MANIFEST_all.yaml` is the union of the per-set manifests (`MANIFEST.yaml`,
`MANIFEST_{pdb,distillation,teddymer,afdb_multimer,msa,cifcore*}.yaml`); regenerate it after
editing one (`python -c 'from structcooker.manifests import write; write()'`, enforced by a
test). Each set's manifest still builds that set alone. What each stage builds is in
[stages](stages.md).

Sizes and times on this cluster (2026-09/10): OFD long cif ≈ 8 × 3 h; teddymer hmmsearch (8 nodes) ≈ 31 h; AFM hmmsearch (8 nodes) ≈ 6 days.
Planning and verification add to that. Measured on OFD long (2026-10-08..12): scanning raw/
(16.4M files) takes ~80 min, hashing ~3,200 files/s, verification ~0.2 CPU s per MSA record
(~8 h for 16M on a 112-core verification job). A 1% sample of everything (`e2e-sample prepare
--count ofd_long=160000 ...`) builds in a few hours.

Large `parallel` searches declare `node_count: 8` and are submitted as an 8-task array; each
task refuses to start on a node whose `/dev/shm` cannot hold process semaphores (joblib would
otherwise run serially there). List such nodes in `DATACOOKER_EXCLUDE_NODES` before
`build-all` so the scheduler does not place tasks on them; scripts/cluster/cssb.env has this
cluster's.

## 5. Check and document the release

```bash
structcooker release-check     # SLURM job -> $BIOMOL_ROOT/release_check.json
structcooker release-docs      # writes $BIOMOL_ROOT/docs/
```

`release-check` checks the tree is complete and self-contained: every stage's output at its
configured path, a shard collection's shards inside its own directory, entry counts agreeing
with `.meta.json` and `.index.tsv`, and a sample of 200 records per DB decoding and passing the
build's checks. `release-docs` writes the release's own README, provenance and licences
(copied from `release_docs/`), the database list, the what-to-copy scenarios with this tree's
sizes (`release_docs/scenarios.yaml`) and a summary of `release_check.json`. Both default to
`db/MANIFEST_all.yaml`; rerun them after a build, a patch or a move.

## 6. Sampled end-to-end check

A full build is days and terabytes. To check that the release DAG still builds end to
end after a change, run every stage on a sample of a few hundred inputs instead
(`structcooker.sample_e2e`). All four commands are login-node safe: `prepare` and `verify`
submit SLURM jobs, and `run` submits stage jobs the way `build-all` does.

```bash
R=/data/shared/cssb_data/e2e_sample/$(date +%Y%m%d)      # any new, empty directory
structcooker e2e-sample prepare $R      # job: pick the sample, build $R as a small BIOMOL_ROOT
structcooker e2e-sample run $R          # build-all on the sample, run state in $R/run
structcooker e2e-sample status $R       # stage phases; rerun `run` to resume after a fix
structcooker e2e-sample verify $R       # job: compare with the release -> $R/verify.json
```

**Sample.** By default 30 PDB entries and 20 entries each from OFD short, long, rna and
disordered (plus 20 disordered template chains), teddymer, and AFM homo- and heterodimers
(`--count pdb=50` etc.). The
selection is a seeded hash order, so the same production data gives the same sample.
On top of that random selection, the sample also includes these cases on purpose:

- the accepted PDB failures (`pdb/cif` in `db/release-exceptions.json`)
- two manually fixed mmCIFs
- two entries from each published PDB-derived DB (train, valid1/2, mpnn, chain), so every filter has something to keep
- the six repaired AFM heterodimers
- the homodimer whose MSA we generated
- two recovered OFD long entries

Their files are linked one by one into `$R/materials` in the release layout: `$R` is itself a
small `BIOMOL_ROOT`.
The reference fastas, disordered work lists, `chain_msa.tsv` and `overrides.json` are
rewritten for the sample. The following are linked whole:

- the CCD
- the SabDab summary
- the seq_id seed
- the reference CCD and chain-selection DBs

No search database is needed: every MSA the build reads is a download or a supplied input.

**Run.** `BIOMOL_ROOT` points at the sample, reader overrides are cleared, and `run`
refuses to start if any config still resolves outside `$R` (or `SEQ_DB_ROOT`). Parallel
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

## 7. Patching one DB

To fix some records of a published DB — built wrong or missing — without rebuilding it:

```bash
structcooker patch distillation/long_cif --files bad_inputs.txt [--nodes 8]   # build op: input files
structcooker patch pdb/cif_attached --keys bad_keys.txt                        # rebuild op: keys
```

The records are built into a separate LMDB next to the target
(`.<name>.build/patch_<stamp>/`), then applied: a single LMDB takes the new values in place; a
shard collection drops those keys from its shards and gains the patch as one more shard. The
index is rewritten from its old rows plus the patch. `--dry-run` writes the scripts without
submitting. A patch is outside the release run's bookkeeping, so rerun `release-check` and
`release-docs` afterwards.
