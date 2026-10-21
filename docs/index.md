# StructCooker

StructCooker builds the **BioMol** database set — structures, MSAs, templates, training and
validation sets as LMDBs, plus the sequence-id and cluster tables — from downloaded sources:
the PDB, the CCD, the OpenFold3 distillation sets, teddymer and the AFDB multimer release.
Every database is one declarative config, `db/<group>/<name>.yaml`; the whole release is the
94-stage DAG in `db/MANIFEST_all.yaml`, built on SLURM by one command.

StructCooker is the biomolecular layer on top of
[DataCooker](https://CSSB-SNU.github.io/DataCooker/), which provides the generic engine:
recipes, the planning-first LMDB pipeline (plan → shard array → merge → index) and the SLURM
executor. How recipes and LMDB builds work is documented there, once; this site covers the
data sets, the configs and how to build and check a release.

## Install

StructCooker vendors DataCooker as a git submodule (`libs/datacooker`) and uses
[pixi](https://pixi.sh):

```bash
git clone --recursive git@github.com:SanggeunParrk/StructCooker.git
cd StructCooker                      # already cloned: git submodule update --init --recursive
pixi install                         # StructCooker + DataCooker (editable), biomol, hmmer, kalign, mmseqs2, ...
pixi run python -c "import structcooker, datacooker, biomol; print('ok')"
```

## Quick start

Everything hangs off one root, `BIOMOL_ROOT`; inputs must already be in place under
`$BIOMOL_ROOT/materials/` (downloads, staging and supplied reference inputs:
[end-to-end build](e2e-build.md) §1–3).

```bash
source scripts/cluster/cssb.env       # BIOMOL_ROOT and the SLURM settings; elsewhere, copy and edit it
pixi run structcooker inspect --strict --manifest db/MANIFEST_all.yaml
pixi run structcooker build-all --manifest db/MANIFEST_all.yaml --workdir /path/to/run   # rerun to resume
pixi run structcooker release-check   # after summary.json says complete: true
pixi run structcooker release-docs
```

Every command submits SLURM jobs and returns; nothing heavy runs where you type it.

## Commands

| Command | Does |
|---|---|
| `structcooker list` | every `db/` config with its op, schema and target |
| `structcooker inspect [--strict] [--manifest M] [NAME]` | preflight: external inputs present, manifest dependencies closed |
| `structcooker download {ccd,sabdab,mmcif,openfold}` | fetch a raw input into `materials/raw` (`mmcif` needs `--yes`; `openfold` prints instructions) |
| `structcooker build NAME [--show/--dry-run]` | build one stage on its own (op inferred from the config) |
| `structcooker build-all --manifest M --workdir RUN` | the release: every stage as verified SLURM jobs, resumable |
| `structcooker patch NAME --keys/--files LIST` | rebuild some records of a published DB and swap them in |
| `structcooker release-check` | check a release tree is complete and self-contained → `release_check.json` |
| `structcooker release-docs` | write the release's own docs into `$BIOMOL_ROOT/docs` (Markdown pages + an `index.html` overview) |
| `structcooker e2e-sample {prepare,run,status,verify} ROOT` | the release DAG on a sample of a few hundred inputs, compared with the release |
| `structcooker foldbench {keys,build,verify}` | FoldBench inference DBs cut from the release by key |
| `structcooker interface {annotate,build,finalize,split,validate}` | the PDB interface-unit DB (outside the `build-all` DAG; see [Interface DB](interface-db.md)) |

## Pages

- [End-to-end build](e2e-build.md) — the runbook: download, stage, reference inputs,
  `build-all`, release check and docs, the sampled check, patching.
- [Stages](stages.md) — what each data set's stages do, their recipes and parameters.
- [Raw materials](raw-materials.md) — the `materials/{raw,reference,intermediate}` tree.
- [seq_id and cluster scheme](seq-id-and-cluster-scheme.md) — one id space, per-DB clusters.
- [FoldBench DBs](foldbench.md) — key-list subsets for FoldBench inference.
- [Interface DB](interface-db.md) — the PDB cut into units of touching molecules (2..N, waters included), with a date split.
- [Development](development.md) — tests, cluster settings, code layout, LMDB storage.
- History — dated reports and superseded guides, kept as records.
