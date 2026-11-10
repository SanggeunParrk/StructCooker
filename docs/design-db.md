# Design DBs

Training sets for backbone design models (La-Proteina style): single-chain OpenFold3 distillation structures, kept only
where the model is confident and the fold is not mostly coil, with the foldseek cluster of every item so a sampler can
balance by structural cluster. There is one DB per OFD monomer set:

| DB | source | stages | 2026-10 hand build |
|---|---|---|---|
| `lmdb/design/db_monomer_long/design_monomer_long.lmdb` | `distillation/long_cif_attached` | `design/monomer_long_*` | 7,969,067 items |
| `lmdb/design/db_monomer_short/design_monomer_short.lmdb` | `distillation/short_cif_attached` | `design/monomer_short_*` | 225,515 items (of 430,418) |

Both are ordinary release DBs: one config per stage under `db_configs/design/`, recipes in `src/structcooker/workflows/design/`,
the functions in `src/structcooker/instructions/design/`, a manifest `db_configs/MANIFEST_design.yaml` that
`MANIFEST_all.yaml` includes. `structcooker build-all` builds them with everything else, and `structcooker e2e-sample`
runs them on the sample (`prepare` needs nothing extra: the stages read the sample's own `distillation/*_cif_attached`).

The earlier scripts in `scripts/design/` (statistics, `cluster_*.sbatch`, `summarize.py`, `build_db.py`, `rebuild_40*.sh`) are
gone; the stages below do the same with the same defaults and give the same statistics (checked on three items of the
2026-10 short DB). Not in the DAG: the PDB complex design DB (dropped, foldseek multimer clustering ran out of memory) and the
Teddymer complexes (already a set of cluster representatives, `teddymer/cif_attached`).

## Stages

For a set `<s>` = `monomer_long` or `monomer_short` (`$R` is `BIOMOL_ROOT`; paths under `materials/intermediate/design/<s>/`
are scratch, rebuildable from the source):

| Stage | Op | Reads | Writes |
|---|---|---|---|
| `design/<s>_stats` | parallel (8-node array) | source attached DB | `stats/part_<n>.tsv` (id, length, mean_plddt, helix, strand, coil, max_coil_run), `backbones/part_<n>.tar` (one N/CA/C/O PDB per item), `stats/part_<n>.failed.tsv` for unreadable records |
| `design/<s>_clusters` | materialize (full node) | `backbones/` | `clusters.tsv`: `representative`, `member` per item |
| `design/<s>_items` | materialize | `stats/`, `clusters.tsv`, source (index) | `$R/metadata/design_<s>_items.tsv`: id, length, mean_plddt, helix, strand, coil, max_coil_run, cluster, cluster_size, pass |
| `design/<s>` | rebuild, schema B | source attached DB, items table | the DB |

- **stats.** The source's sorted keys are cut into parts of `chunk_keys` (5000); a part is one work item. Per item the
  first assembly is read as a `CIFMolAttached`: `length` is the number of standard amino-acid residues with N, CA, C and O,
  `mean_plddt` the mean B factor of all atoms (OpenFold stores pLDDT there), the secondary-structure fractions and the longest
  coil run come from the vectorised DSSP of `instructions/design/ss.py` (checked against `mkdssp` by
  `scripts/design/check_ss.py`). The part's `.tsv` is written last, so a resumed run lists only parts without one.
- **clusters.** `foldseek createdb` on all backbone archives, `foldseek cluster` with its defaults and `-c 0.9`, `createtsv`.
  `n_jobs` (112) is the thread count and makes the job a full node. An item without any backbone residue has no structure to
  cluster; the items stage makes it its own cluster.
- **items.** Joins the three and sets `pass`. It refuses to build if a source key has no statistics, if more than `max_failed`
  records failed in stats (0), or if an item with a backbone is in no cluster. `cluster_size` counts the cluster's members
  among all items, passing or not. Thresholds only change this stage and the DB, not the foldseek run.
- **DB.** Rebuilds the source keeping the items with `pass = 1`, with the same keys and the same
  `{assembly key: {"cifmol_attached_dict": ...}}` layout, adding to each record's CIFMolAttached `metadata`:

  | key | value |
  |---|---|
  | `source` | `distillation_long` / `distillation_short` (`parameters.source`) |
  | `length` | residues with a complete backbone |
  | `mean_plddt` | mean B factor |
  | `ss_helix`, `ss_strand`, `ss_coil` | fractions of the assessed residues |
  | `ss_max_coil_run` | longest run of consecutive coil residues |
  | `fs_cluster` | id (source key) of the foldseek representative |
  | `fs_cluster_size` | members of that cluster in the source |

  The items table travels to the workers as a node-local LMDB (as the sequence metadata of the attach stages does).

## Parameters (all in the configs)

| Parameter | Value | Config |
|---|---|---|
| `min_plddt` | 80.0 | `<s>_items` |
| `min_length`, `max_length` | 32, 512 | `<s>_items` |
| `max_coil` | 0.5 (coil fraction strictly below) | `<s>_items` |
| `max_coil_run` | 20 | `<s>_items` |
| `apply_filters` | true (false: every item passes, for a set filtered when it was made) | `<s>_items` |
| `max_failed` | 0 | `<s>_items` |
| `coverage` | 0.9 (foldseek `-c`) | `<s>_clusters` |
| `split_memory_limit` | 400G (foldseek) | `<s>_clusters` |
| `chunk_keys` | 5000 | `<s>_stats` |
| `parameters.source` | `distillation_<long|short>` | `<s>` |

## Build

```bash
structcooker build-all --manifest db_configs/MANIFEST_design.yaml --workdir <run dir>   # this set, with the roots it needs
structcooker build design/monomer_short_items                                          # one stage on its own
structcooker build design/monomer_long --show                                          # resolved invocation
```

`foldseek` must be on the jobs' PATH (it is in the pixi environment; `structcooker inspect` lists the tools).

Sizes to plan for (long): the backbone archives are plain PDB text, about 80 KB per item, so several hundred GB under
`materials/intermediate/design/monomer_long/backbones/`; every stage that reads them fingerprints them, and the foldseek
tmp directory is several times larger while `clusters` runs (it is removed afterwards). The long clustering needs a full CPU
node (the nodes have about 503 GB; `split_memory_limit` is 400G). The `items` stage fingerprints the source DB (it checks
the statistics against the source's index), as every stage that reads a DB does.

## In the release DAG

- `MANIFEST_design.yaml` lists the stages and, as roots, the distillation stages they read (with the dependencies those
  have in `MANIFEST_distillation.yaml`); `MANIFEST_all.yaml` is regenerated from it
  (`python -c 'from structcooker.release.manifests import write; write()'`, a test enforces it).
- **Sampled end-to-end.** `e2e-sample run` builds the eight stages on the sample. `verify` compares every LMDB with the release
  record of the same key; a DB without a release counterpart is reported as `"release": "absent"` and is not a failure (schema
  B is not in `EXACT_SCHEMAS`), the projections are listed by file count. Nothing needs a release-exceptions entry: those
  name inputs a stage may fail on, and the design stages fail on none (`max_failed: 0`). With the default 20 short and 20 long
  items the DB is small but not empty (about half of the short items pass).
- **release-check** requires `.meta.json` and `.index.tsv` to agree with the DB, which the rebuild writes.

## Developer tools

`scripts/design/check_ss.py` (DSSP against `mkdssp`), `probe_monomers.py` (small end-to-end trial on the short set with
`foldseek easy-cluster`), `validate_db.py DB SOURCE` (a built DB against its source and fresh statistics), and
`run_py.sbatch` to run any of them as a job. They are not pipeline stages.
