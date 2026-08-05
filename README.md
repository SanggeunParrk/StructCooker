# StructCooker

StructCooker reproduces the **BioMol** LMDB database set — the CIF / MSA / template /
training DBs used to train structure-prediction models — from raw inputs, using
**declarative configs** on the [DataCooker](https://CSSB-SNU.github.io/DataCooker/)
engine. Every database is one config; you build it with one command.

> **One artifact = one `db/**/<name>.yaml`.** `structcooker build <name>` hands that
> config to the DataCooker pipeline, which shards the work across nodes, runs the
> recipe, merges, and indexes — with **Ray** owning memory (no hand-tuned
> `n_jobs`/`mem`/`shards`). Nothing about a database lives outside its config.

## Quickstart

```bash
# 1. clone with the engine (vendored as a submodule)
git clone --recursive <repo-url> && cd StructCooker
# already cloned? git submodule update --init --recursive

# 2. environment (pixi installs StructCooker + DataCooker + biomol + hmmer/kalign/…)
pixi install && pixi shell

# 3. point the two roots at your machine (portable across servers)
export DATA_ROOT=/path/to/raw/downloads   # mmCIF, CCD, OpenFold, … (read)
export OUTPUT_ROOT=/path/to/reproduced/db # where built LMDBs go   (write)

# 4. (optional) to reproduce the EXISTING BioMol set identically, seed the seq_id map
#    -- seq_id is an assigned counter, so from scratch you'd get a different (but
#    self-consistent) id space. Omit this to build a fresh set.
python -c "from huggingface_hub import hf_hub_download; \
  hf_hub_download('biomol/seq-id-map','seq_id_map.tsv.gz',repo_type='dataset',local_dir='$DATA_ROOT/metadata')"
gunzip -f $DATA_ROOT/metadata/seq_id_map.tsv.gz
export SEQID_SEED=$DATA_ROOT/metadata/seq_id_map.tsv

# 5. fetch the raw external inputs (seq_id_map is the HF seed above, not here)
structcooker download ccd            # wwPDB CCD -> OUTPUT_ROOT/materials/raw/ccd
structcooker download sabdab         # SabDab antibody summary (seq_cluster input)
structcooker download mmcif --yes    # full wwPDB mmCIF (~90 GB+) -- needs --yes
structcooker download openfold       # TB-scale: prints portal instructions, no auto-fetch

# 6. preflight: is everything a build needs already in place?
structcooker inspect                 # per-node READY/BLOCKED + missing external inputs

# 7. see what you can build
structcooker list

# 8a. build one database (op auto-inferred from the config)
structcooker build pdb/cif           # raw mmCIF   -> pdb/cif        (build)
structcooker build metadata/seq_id_map  # cif fasta -> seq_id map    (materialize)
structcooker build pdb/cif_attached  # + metadata  -> pdb/cif_attached (rebuild)

# 8b. or reproduce the whole DAG in dependency order, incrementally
structcooker build-all               # skips whatever is already built
```

Run `structcooker inspect` first: it resolves every build-all node's inputs, checks the
external tools + the two env roots, and lists exactly which raw inputs are missing (each
tagged with the `structcooker download` target that fetches it, or "provide" for a
tool/lab output) — so you know a build-all can actually complete before you launch it.

`structcooker build` auto-infers the op from each config: **build/rebuild** run the
planning-first SLURM pipeline (tiers → merge → index, afterok-chained); **materialize/
extract** (metadata projections that write a TSV/fasta, not an LMDB) run a single
workflow job. `build-all` wires every config into one DAG (`db/MANIFEST.yaml`) and
submits them SLURM-ordered, skipping already-built outputs. You need a SLURM cluster
and the raw inputs on disk. See [docs/roadmap.md](docs/roadmap.md) for raw inputs.

**Reproducing production identically vs. a fresh set.** Most DBs rebuild decode-level
identical to production from raw inputs. The exception is `seq_id_map`: seq_id is an
assigned running counter, not a hash, so an unseeded build mints a *different* (but
internally coherent) id space, which then threads through everything keyed by seq_id.
Set `SEQID_SEED` to the published map (above) to match production; leave it unset for a
fresh set. `seq_cluster` matches production given the same corpus (`SEQCLUSTER_FASTA`)
+ mmseqs2 version — it is deterministic, so it needs no seed.

## What you can build

`db/STATUS.md` is the ledger: which databases have a clean config today (✅) versus
still need porting (🔴) or are deferred (⏸️). [docs/roadmap.md](docs/roadmap.md) has
the full dependency DAG for the three deliverables — the PDB set + training DBs,
the OpenFold distillation sets, and the MPNN training view.

Every config resolves its paths from two env vars — `DATA_ROOT` (raw inputs, read)
and `OUTPUT_ROOT` (reproduced DBs, write) — so the same configs run on any server by
pointing those two at the right places (they default to this cluster's layout).
Correctness bar: a rebuilt DB is decode-level identical to production.

## How it fits together

```
db/**/*.yaml                     the build interface — one config per database
src/structcooker/
  cli.py                         structcooker build | list
  schemas.py                     value-type catalogue (schema tag on each config)
  instructions/{readers,transforms,writers}/   reusable transform primitives
  workflows/{ingest,filters,metadata,exports,search,analysis}/   recipe modules
libs/datacooker/                 the engine (planning-first, Ray fan-out)
docs/                            roadmap, per-topic references
```

A config names a `recipe` and its `reader`/`writer` hooks (all in
`src/structcooker/…`) plus a source and target; `structcooker build` strips the
`schema:` tag and hands the rest to the engine.

## Docs

- [Getting started](docs_src/getting-started.md) — clone, install, smoke test
- [docs/roadmap.md](docs/roadmap.md) — the reproduction plan + DAG
- [db/STATUS.md](db/STATUS.md) — per-database status ledger
- [DataCooker](https://CSSB-SNU.github.io/DataCooker/) — the engine
