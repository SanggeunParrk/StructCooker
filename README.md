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

# 4. see what you can build
structcooker list

# 4. build one database (submits the planning-first SLURM pipeline)
structcooker build pdb/cif          # raw mmCIF  -> pdb/cif
structcooker build pdb/cif_attached # + metadata -> pdb/cif_attached
structcooker build msa/a3m_d16k     # a3m depth-capped to 16000
```

`build` submits a SLURM job array (tiers → merge → index, afterok-chained), so you
need a SLURM cluster and the raw inputs on disk. See
[docs/roadmap.md](docs/roadmap.md) for how raw inputs are obtained.

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
