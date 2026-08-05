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
| `chain/cif_chain` | rebuild | C | cif_pdb_attached `.index.tsv` (split→explode per chain; verified 333/333 decode-identical) |

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

## 🟢 Metadata + projection ops (new op class)

`structcooker build` now auto-infers **materialize** (recipe → TSV/fasta) and
**extract** (LMDB → TSV) ops alongside build/rebuild, routed to
`datacooker.cli.workflow` as a single afterok-chainable job. The PDB metadata chain
is ported under `db/metadata/`:

| config | op | notes |
|---|---|---|
| `cif_fasta` / `cif_metadata` | extract | cif_pdb → fasta / metadata TSV (recipes verified on prod records) |
| `seq_id_map` | materialize | **optional seed exception** — seq_id is an assigned counter, so from-scratch ≠ production; seed with the published map (HF `biomol/seq-id-map`, `SEQID_SEED`) to match, else fresh ids. Validated. |
| `seq_cluster40` / `seq_cluster30` | materialize | mmseqs2 (antibodies via cd-hit); deterministic given corpus+params+version. Corpus = `SEQCLUSTER_FASTA` (prod used the pdb+distillation union); SabDab is an external input |
| `interacting_seq_ids` / `interacting_seq_clusters` | extract / materialize | interface partners for valid-2 dedup |

### Key-list sizing gap — CLOSED
A build with `keyed: true` sizes items by **count** (uniform), not input `st_size`,
so key-list builds (items are seqids/chain-ids, resolved by lookup) plan cleanly.

### Template hmm pipeline — ABSORBED (4th op: `parallel`)
The full template pipeline is now on the clean surface, wired in MANIFEST order
`hmmsearch → template_candidates → template_seqs → seqid_template_mols → chain_template`:

| config | op | notes |
|---|---|---|
| `template/hmmsearch` | parallel | Phase 0 — hmmbuild(query a3m→HMM) + hmmsearch(vs L-chain DB). `split_recipe` → `datacooker.cli.workflow parallel-run`. Idempotent |
| `template/template_candidates` | materialize | Phase 1/2 — seqid_to_chains, chain_to_templates, reduced_hmm (verbatim port; logic verified on synthetic inputs) |
| `template/template_seqs` | materialize | Phase 3 seq maps (seqid_to_seq, chain_to_seq) — supersedes the old cif_chain_seq LMDB |
| `template/seqid_template_mols` / `chain_template` | build (keyed) | Phase 3/4, schema D |

`msa/msa_rna`, `valid/valid1`, `valid/valid1_attach`, `valid/valid2` are also ported.

### Minor derivations — DONE
The small list/fasta projections feeding the hmm pipeline are ported as materialize/
parallel ops: `metadata/pdb_polypeptide_L` (L-chain fasta = hmmsearch template DB),
`metadata/template_chain_filelist` (protein-chain work list), `metadata/seqid_template_filelist`
(Phase 3 keyed-build item list), `template/msa_wo_lower` (a3m insertion-strip). The full
**31-node MANIFEST** topo-sorts end to end.

## ✅ PDB reproduction — complete on the clean surface

Every step from raw mmCIF to the training/validation DBs is a `db/**` config; nothing
runs off a legacy script. What is *not* reproduced here (by design, provided like the
raw mmCIF/CCD downloads): the external-tool / raw inputs — **SabDab** (antibody summary),
**SignalP** outputs, and the raw per-chain **a3m** MSAs. The legacy `rekey_seq_id_db.py`
is a one-off width migration (`P0000007` → `P…020d`) for *pre-existing* DBs; a
from-scratch build already emits 20-width ids (`_SEQ_ID_WIDTH`), so it is not a
reproduction step.

### Manual mmCIF fixes — a required pre-ingest step (provided input)

53 PDB entries (`db/pdb/manual_cif_fixes.txt`) error out / build wrongly from the current
wwPDB mmCIF — mostly NMR ensembles whose non-polymer ligand is re-numbered per model,
breaking the atom→scheme match (verified on `1ai0`: `IPH` at `auth_seq` 22 in some models,
31 in others, vs a single scheme number). The production build substituted an older
known-good cif for each before ingest (legacy `scripts/manually_fix_cif.py`, source
`BioMolDB_2024Oct21`). Ported as `structcooker fix-cif --source <corrected-cif-dir>`; run
it before `build pdb/cif`, and `inspect` reports the state. **Not** a port regression (the
cif logic is byte-identical across the whole repo history) and **not** a CCD difference.
Full write-up: [docs/manual-cif-fixes.md](../docs/manual-cif-fixes.md).

**Current decision:** the corrected snapshot is not on this cluster, so the accepted
`cif_pdb.lmdb` (233,579 entries) is built **without** the substitution — 25 of the 53 are
absent (errored), 28 build from the current mmCIF (present but not production-substituted).
`fix-cif` is ready to fold all 53 in production-faithfully once the snapshot is obtained.

The OpenFold3 distillation sets. **Recipes exist** (`workflows/ingest/openfold_*`);
configs are ported onto the clean `db/distillation/` surface (env-var paths, schema tag,
ops knobs dropped — same shape as long_cif). Production stores the *downstream* DBs
(`cif_*_attached`, capped `msa_*_{d2k,d512,d5k}`, `template_topn_*`/`_disordered`), so the
full chain is base → derived.

**✅ ported + resolve + small-validated (21 configs):**

| kind | schema | configs |
|---|---|---|
| structure (base cif) | A | `{long,short,rna,disordered}_cif` |
| msa (base, depth 65536) | E | `{long,short,rna,disordered}_msa` |
| msa depth-caps | E | `{long,disordered}_msa_{d2k,d512}`, `{short,rna}_msa_{d2k,d512,d5k}` |
| template (top-N) | D | `{long,short}_template_topn` |
| template (disordered) | D | `disordered_template` |

Base-structure pipeline validated in-process (Ray-free): one real `structure.npz`
(`MGYP003648360693`) → recipe → schema-A valid → serialize/deserialize round-trips.

**✅ `cif_*_attached` (schema B) — ported (25 configs total).** Previously blocked on a
rewrap step + missing metadata; both resolved. The base `cif_{set}` already emits the
wrapped `assembly_dict`/`metadata_dict` layout (the openfold_structure recipe), so no
separate rewrap is needed, and `seq_id_map` / `seq_cluster40` are now on the clean
surface (see the metadata section). `cif_{long,short,rna,disordered}_attached` are plain
A→B rebuilds of the base + the shared attach recipe (proven by pdb/cif_attached).

By project rule these configs are code-complete only: validate on a small sample vs the
existing production DBs, **do not full-build** (they already exist).
