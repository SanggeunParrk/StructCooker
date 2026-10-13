# StructCooker

StructCooker builds the **BioMol** database set — the structure, MSA, template, training and
validation LMDBs used to train structure-prediction models — from downloaded sources (PDB,
CCD, OpenFold3 distillation, teddymer, AFDB multimer), on the
[DataCooker](https://CSSB-SNU.github.io/DataCooker/) engine. Every database is one config,
`db/<group>/<name>.yaml`; the release is the 94-stage DAG in `db/MANIFEST_all.yaml`, built on
SLURM by `structcooker build-all`.

```bash
git clone --recursive <repo-url> && cd StructCooker && pixi install
source scripts/cluster/cssb.env       # BIOMOL_ROOT and SLURM settings (copy it for another cluster)
pixi run structcooker inspect --strict --manifest db/MANIFEST_all.yaml
pixi run structcooker build-all --manifest db/MANIFEST_all.yaml --workdir /path/to/run
```

- [docs/index.md](docs/index.md) — what it is, install, quick start, every command.
- [docs/e2e-build.md](docs/e2e-build.md) — the build runbook: downloads, staging, reference
  inputs, `build-all`, release check, patching.
- [docs/stages.md](docs/stages.md) — what each stage builds and with which parameters.
- [db/STATUS.md](db/STATUS.md) — what the current release holds, its known gaps and checks.

Sequence ids are assigned counters: keep `materials/reference/seq_id_map_seed.tsv` (and the
other supplied reference inputs, [docs/e2e-build.md §3](docs/e2e-build.md#3-reference-inputs))
to reproduce this release's ids; without them a build is self-consistent but different.

Code: MIT. The data keep their sources' terms (`release_docs/licenses.md`).
