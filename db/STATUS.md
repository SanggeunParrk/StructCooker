# db/ migration ledger

Every non-distillation BioMol DB and the status of its declarative `db/*.yaml`.
`structcooker build <name>` works only for ✅ rows; 🔴 rows need engine work first
(a build whose items are **keys** can't be sized by input `st_size`, and script-only
builds must be absorbed into a datacooker recipe/op before they get a config).

## ✅ Clean — config written, planning-first works

| db config | op | schema | source of truth |
|---|---|---|---|
| `pdb/cif` | build | A | `mmcif_files_latest` (real .cif files → st_size) |
| `pdb/cif_attached` | rebuild | B | cif_pdb `.index.tsv` |
| `train/train_20210930` | rebuild | B | cif_pdb_attached `.index.tsv` |
| `train/train_20260301` | rebuild | B | cif_pdb_attached `.index.tsv` |
| `train/train_20210930_mpnn` | rebuild | H | train_20210930 `.index.tsv` |
| `msa/a3m` | build | E | real `*.a3m` files → st_size |
| `msa/a3m_d16k` / `a3m_d2k` / `a3m_d512` | rebuild | E | a3m `.index.tsv` (depth cap 16000/2000/512) |
| `ccd/ccd` | build | G | component `.cif` files → st_size |

Replaces the 10 hand-tuned `cif_pdb_{light,medium,large,xlarge,huge,monster,smoke,…}`
variants with `pdb/cif.yaml` alone. Schema **G** (CCD chem-component) added for `ccd`.

Schema **H** (pruned CIFMol, one key per assembly-model-altloc) added for the first
*model-input* DB, `train/train_20210930_mpnn` — a view of a CIF DB carrying only the
features one model reads. It is the one rebuild that is 1:N rather than 1:1
(`explode_entries`), so a key is a training sample. See `docs/mpnn_pruning.md` for the
size measurements and for why its `index_table` must be read through
`instructions.transforms.pruning.load_pruned_cifmol`.

### MSA depth-cap bypass absorbed (was a 🔴 row)
`scripts/maintenance/lightweight_msa.py` + `build_msa_long_d2k.sh` are replaced by the
`workflows/exports/cap_msa_depth.py` rebuild recipe (+ `adapt_msa_for_cap` adapter):
`db/msa/a3m_{d16k,d2k,d512}.yaml` differ only in `parameters.max_depth`. Verified on
a3m_d16k (178,249 keys, key-set complete, deep records capped to 16000, standard-frame
output). The recipe caps `msa_dict` in place; the adapter renames the input so the
same-named output does not silently return the *uncapped* record.

### Ray is the build/rebuild fan-out backend — Ray-only (libs/datacooker submodule)
The **backpressure engine open item is closed**: `build_lmdb`/`rebuild_lmdb` no longer
predict memory (`n_jobs × E × max_bytes`) and hand-tune per DB — they fan out over **Ray**,
whose node memory monitor throttles+auto-retries tasks by *measured* memory (admission
control). There is **no backend choice**: the `parallel_backend` config knob is gone, the
fork (`multiprocessing`) / loky branches are deleted, and `parallel_backend` /
`maxtasksperchild` / `task_timeout` are accepted-but-ignored for config back-compat only.
See `datacooker/_ray.py` and [[datacooker-memory-aware-engine]].

* **Measured, not predicted.** Benchmark on the 22,281 deepest a3m records (the slice that
  forced fork down to `n_jobs=16`): fork **7 rec/s @ 53 G** (memory-capped) vs shipped Ray
  **74 rec/s @ 225 G** (actor prototype hit 94 rec/s) — **~10× faster**, no tuning, and
  output **deep-equal identical** (2229/2229 sample + 300/300 no-knob run). Ray never had
  to throttle: 225 G is 46 % of the 491 G node.
* **Recipe philosophy intact.** Ray runs the *same* `_rebuild_worker → rebuild_entry →
  parse_dict`; the recipe is backend-blind (that decoupling is *why* the swap was a
  one-liner and why output is identical). Only the fan-out layer changed.
* **State shipped once via plasma** (`ray.put`) → each worker reads zero-copy from
  `/dev/shm`, replicating the fork path's copy-on-write benefit without predicting memory.
* **Native Ray runtime policy** — no custom `/dev/shm` scoping or rmtree: Ray isolates each
  `ray.init()` in its own `session_*` dir (concurrent SLURM shards don't collide) and cleans
  it on graceful shutdown. A SIGKILL leak is bounded (per-session plasma file) and left to
  ops (`ray stop` / epilog), not the library.
* **Scope:** Ray governs the **within-node** fan-out; cross-node tier/shard splitting still
  applies (planner unchanged, its `n_jobs` now advisory). `iter_parallel_chunks` (joblib)
  is **kept** — `extract.py` (edge_node, with its SIGALRM per-record watchdog) and
  `processing/batch.py` still use it; those are separate paths, not the OOM-prone one.
* Still carry forward (independent of backend): **`--chunk-size` bounded** (parent drain
  buffer), **MSA schema E 40→100** (dense-cap sizing feeds tier splitting), **extract
  SIGALRM watchdog**, **executor `cpu-long-q`** (only uncapped cpu QOS).

**Dead-code cleanup done (Ray-only is the single path).** Removed everywhere: the fork
branch (`gc.freeze`/`gc.disable`/`maxtasksperchild`) in `_parallel.py`; `maxtasksperchild`
+ `task_timeout` from `build_lmdb`/`rebuild_lmdb`/`_parallel_write`; `parallel_backend`
knob (params + all config lines); the `Tier.maxtasksperchild` field + planner recycle
logic; `--max-tasks-per-child`/`--task-timeout` on the build/rebuild CLI; `cfg.task_timeout`
in the pipeline. `task_timeout` survives only on the separate `extract` path (its SIGALRM
watchdog). Full `libs/datacooker/tests` suite **53/53 green** under ray-only.

> **Ops gotcha:** build/rebuild now spin a Ray cluster, so their tests (and any real build)
> must run on a **compute node** — never the login node. Run the suite with
> `PYTHONPATH=libs/datacooker` so Ray workers (fresh processes) can import test-module
> recipes that the fork backend used to inherit via COW.

### Superseded old configs deleted (git-recoverable)
`configs/ingest/cif_pdb_{light,medium,large,large_fast,xlarge,xlarge_fast,huge,monster}`,
`cif_lmdb`, `cif_lmdb_refactor_test`, `filters/train_filter{,_20260224}`, `ingest/a3m_lmdb`,
`ingest/ccd_lmdb` (14). Untracked-and-superseded, left for now (no git safety net):
`cif_pdb_large_retry`, `cif_pdb_smoke`, `metadata/attach_cifpdb_revisit`, `ingest/ccd_lmdb_20260718`.
`logs/` cleared (5894 files / 393 MB).

## 🟡 Multi-stage — belongs to the step-4 pipeline DAG, not a single config

| db | why |
|---|---|
| `valid/valid1` (stage 1) → attach → `valid/valid2` (stage 2) | stage 2 dedups sequences vs **train/valid_1 fasta + seq clusters** — cross-set deps (`validation_stage2*`), needs the orchestration DAG, not a lone rebuild config |

## 🔴 Needs engine work (can't just write YAML)

| db | current builder | why not clean |
|---|---|---|
| `chain/cif_chain` | `scripts/maintenance/build_cif_chain.py` | script, not a datacooker config — absorb into a recipe |
| `chain/cif_chain_seq` | `scripts/maintenance/build_chain_seq.py` | ditto |
| `template/seqid_mols` (Phase 3) | `ingest/seqid_template_mols` | **items are seqids, not files** → planner st_size sizing N/A; also multi-input (metadata + cif_chain lookup) |
| `template/chain` (Phase 4) | `ingest/chain_template_from_seqid` | items are chain keys, not files — same sizing gap |
| `template/pdb` (lower) | `scripts/maintenance/rekey_*` | uppercase→lowercase rekey — absorb as a datacooker `rekey` op |
| `template/candidates` (Phase 2) | `scripts/maintenance/precompute_template_candidates.py` + `cat` | custom py + shell — the plan's "bypass #1" |
| `msa/msa_rna*` | (rna msa scripts) | same — port to the cap_msa_depth pattern once RNA MSA is in scope |

### Engine gap to close (feeds datacooker step 1)
Planning-first `build` sizes items by input-file `st_size`. Key-list builds
(template Phase 3/4) have no input file — their cost is the **output/lookup** size.
Options: size uniformly (n items × mean), or read the source `cif_chain`/seqid index
for per-key bytes. Until then these stay on their existing submits.

## ⏸️ Deferred — code/config only, build later (per user)

`distillation/{long,short,rna,disordered}`: `cif_*_attached`, `msa_*_{d2k,d512,d5k}`,
`template_topn_*`, `template_disordered`. Write configs + recipes, tag build deferred.
