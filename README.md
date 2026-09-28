# StructCooker

StructCooker reproduces the **BioMol** LMDB database set — the CIF / MSA / template /
training DBs used to train structure-prediction models — from raw inputs, using
**declarative configs** on the [DataCooker](https://CSSB-SNU.github.io/DataCooker/)
engine. Every database is one config; you build it with one command.

> **One artifact = one `db/**/<name>.yaml`.** `structcooker build <name>` hands that
> config to the DataCooker pipeline, which shards the work across nodes, runs the
> recipe, merges, and indexes — with **Ray** respecting allocated CPUs, gradually admitting work and
> pausing submission under measured memory pressure. Nothing about a database lives outside its config.

## Quickstart

LMDB finalization uses a [5M-item threshold](docs/lmdb-storage-policy.md): larger
databases keep their shards and publish a manifest instead of copying into one LMDB.

```bash
# 1. clone with the engine (vendored as a submodule)
git clone --recursive <repo-url> && cd StructCooker
# already cloned? git submodule update --init --recursive

# 2. environment (pixi installs StructCooker + DataCooker + biomol + hmmer/kalign/…)
pixi install && pixi shell

# 3. choose input and output roots (defaults match the existing cluster layout)
export DATA_ROOT=/path/to/inputs
export OUTPUT_ROOT=/path/to/reproduced/db # defaults to $DATA_ROOT/BioMol_clean
# Optional input relocation, shared by download helpers and build configs:
export MMCIF_ROOT=$DATA_ROOT/BioMol/materials/raw/cif
export DISTILLATION_ROOT=$DATA_ROOT/BioMol/materials/openfold_distillation

# 4. provide sequence references, SignalP results, and reproduction references.
#    See docs/data-provenance.md; downloads alone do not supply all inputs.

# 5. download into an input directory you own (these commands write raw files)
structcooker download ccd
structcooker download sabdab
structcooker download mmcif --yes
structcooker download openfold       # prints portal instructions; no auto-fetch

# 6. verify input presence and declared dependencies before building
structcooker inspect --manifest db/MANIFEST_cifcore.yaml --strict

# 7. see what you can build
structcooker list

# Run build controllers in a SLURM compute allocation on this cluster.
# 8a. build one database (op auto-inferred from the config)
structcooker build pdb/cif           # raw mmCIF   -> pdb/cif        (build)
structcooker build metadata/seq_id_map  # cif fasta -> seq_id map    (materialize)
structcooker build pdb/cif_attached  # + metadata  -> pdb/cif_attached (rebuild)

# 8b. PDB CIF + supplied MSA: verified resume, input hashes and failure policy
# See docs/pdb-production.md for the SLURM wrapper and required reference inputs.
structcooker pdb-build --run-dir /path/to/durable/pdb-run

# Legacy whole-DAG entry (existence-based skipping, without verified receipts)
structcooker build-all               # skips whatever is already built
```

`inspect --strict` checks required input paths and the manifest's producer/dependency
relationships; missing inputs or invalid dependencies produce a nonzero exit code.
It also lists available tools. Input presence is not a scientific validity check or
proof that an external tool will run successfully.

`structcooker build` auto-infers the op from each config: **build/rebuild** run the
planning-first SLURM pipeline (shards → merge → index, afterok-chained); **materialize/
extract** (metadata projections that write a TSV/fasta, not an LMDB) run a single
workflow job. `build-all` wires every config into one DAG (`db/MANIFEST.yaml`) and
waits for upstream terminal jobs before planning their dependents, skipping already-built
outputs whose upstream nodes were also skipped. Keep the command running until it exits;
it waits for final jobs as well as upstream jobs. A failed or unconfirmable job causes a nonzero exit. You need a SLURM cluster
and the raw inputs on disk. See [docs/roadmap.md](docs/roadmap.md) for raw inputs.

**Reproduction inputs matter.** Sequence IDs are assigned counters. Provide
`DATA_ROOT/reference/seq_id_map.tsv` to preserve existing IDs; its absence creates
a fresh ID space. The clustering recipes also require the supplied union corpus
`DATA_ROOT/reference/seqcluster_corpus.fasta`. Attachment and validation read the
clusters generated under `OUTPUT_ROOT/metadata/` by default. Set `SEQ_CLUSTER30_PATH`
and `SEQ_CLUSTER40_PATH` to supplied cluster files when reproducing an existing
cluster space. Preserve input versions and parameters as well as those files.

For PDB operations, use the [verified build and recovery guide](docs/pdb-production.md).
The completed CIF recovery has full key-set checks and sampled payload comparisons,
with documented chemistry differences and five rejected source entries. See the
[CIF verification report](docs/recovery-2026-09-10.md) and
[MSA reproduction report](docs/msa-reproduction.md) for the bulk-output evidence and
coverage limitations. A fresh full-scale run through the new verified entry point,
template bulk reproduction and a fresh installation have not been verified.

## What you can build

`db/STATUS.md` is the ledger: which databases have a clean config today (✅) versus
still need porting (🔴) or are deferred (⏸️). [docs/roadmap.md](docs/roadmap.md) has
the full dependency DAG for the three deliverables — the PDB set + training DBs,
the OpenFold distillation sets, and the MPNN training view.

Two predicted-complex sets sit alongside them, each with its own manifest:
**teddymer** (`db/MANIFEST_teddymer.yaml` -- 510,454 TED domain-pair dimers, MSAs by
MMseqs2) and **AFDB multimer** (`db/MANIFEST_afdb_multimer.yaml` -- 1,750,755 homodimers
and 80,248 heterodimers, the release's own MSAs). All DBs share one `seq_id` space and
cluster per DB (`cPDB_` / `cOFD_` / `cTDM_` / `cAFM_`); see
[docs/seq-id-and-cluster-scheme.md](docs/seq-id-and-cluster-scheme.md).

Build configs use `DATA_ROOT` and `OUTPUT_ROOT`, with optional input-location
overrides described in [data provenance](docs/data-provenance.md). Existing shared
input paths remain the defaults. A custom `DATA_ROOT` also relocates the default
output root; it no longer silently writes to the cluster's default output directory.

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
- [docs/data-provenance.md](docs/data-provenance.md) — input locations, snapshot references, and the limits of reproduction
- [docs/raw-materials.md](docs/raw-materials.md) — what each staged raw tree holds and where it came from
- [docs/seq-id-and-cluster-scheme.md](docs/seq-id-and-cluster-scheme.md) — one shared seq_id space, per-DB clustering
- [docs/build-all.md](docs/build-all.md) — the build-all DAG + submit/skip/afterok mechanism, diagrammed
- [docs/manual-cif-fixes.md](docs/manual-cif-fixes.md) — historical substitutions, applied only when reproducing that snapshot
- [docs/roadmap.md](docs/roadmap.md) — the reproduction plan + DAG
- [db/STATUS.md](db/STATUS.md) — per-database status ledger
- [DataCooker](https://CSSB-SNU.github.io/DataCooker/) — the engine

The [MSA workflow](docs/msa-reproduction.md) documents supplied alignment inputs,
RNA snapshot differences, and protein/RNA depth caps.

## Development checks

The [quality follow-up](docs/quality-2026-09-10.md) records the latest checks.
See [development.md](docs/development.md) for Ruff, Pyright, and the full test suite
on a SLURM compute node. Recovery logs and snapshots are evidence, not required
imports or build entry points.

For finite native SLURM release submission, verified reuse and restart recovery,
see [Verified release builds](docs/release-build.md). `build-all` now requires a
persistent `--workdir`; `release-build` defaults to the distillation manifest.
