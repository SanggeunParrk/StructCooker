# Getting Started

StructCooker depends on [DataCooker](https://CSSB-SNU.github.io/DataCooker/),
which is vendored under `libs/datacooker` as a submodule.

## Clone (with the engine)

```bash
git clone --recursive git@github.com:SanggeunParrk/StructCooker.git
cd StructCooker
# or, if already cloned:
git submodule update --init --recursive
```

## Install the environment

StructCooker uses [pixi](https://pixi.sh):

```bash
pixi install
pixi shell
```

This installs StructCooker (editable), the DataCooker engine (editable, from
`libs/datacooker`), `biomol`, and the bioinformatics tools used by the
workflows (hmmer, kalign, cd-hit, anarci, ...).

## Smoke test

```bash
pixi run python -c "import structcooker, datacooker, biomol; print('ok')"
```

## Next steps

- [OpenFold3 distillation ingest](tutorials/openfold-distillation.md)
- New to the engine? Start with the
  [DataCooker concepts](https://CSSB-SNU.github.io/DataCooker/) and
  [Build an LMDB](https://CSSB-SNU.github.io/DataCooker/tutorials/build_lmdb/).
