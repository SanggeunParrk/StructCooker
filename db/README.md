# `db/` — one BioMol DB, one declarative config

Each `db/**/<name>.yaml` fully describes how to build **one** BioMol LMDB. It is a
[datacooker](../libs/datacooker) *engine config* (recipe / reader / writer + source +
target) plus a single **`schema:`** tag. Nothing else — no `n_jobs`, `chunk_size`,
`file_list`, or size tiers. Those are derived, not authored.

```
structcooker build <name>        # plan + build (op auto-inferred), afterok-chained
structcooker build <name> --dry-run   # plan + write sbatch, do not submit
structcooker build <name> --show      # print resolved schema/E + invocation, exit
structcooker list                # every db/*.yaml with op, schema, target
```

## How a build runs

Finalization follows the [5M storage policy](../docs/lmdb-storage-policy.md): at most
5,000,000 output items are merged; larger outputs retain immutable `shard_lmdb`
collections. This decision uses completed output counts, including exploded keys.

For PDB CIF and supplied MSA, use `structcooker pdb-build --run-dir <durable-path>`
with `MANIFEST_pdb.yaml`. It adds content-based reuse, durable SLURM restart state,
failure accounting and mandatory output validation. See the
[operating guide](../docs/pdb-production.md); the individual `build` interface below
does not provide that completion contract.

`structcooker build` reads the `schema:` tag, looks up that schema's measured
**expansion factor E** (peak-live-mem ÷ stored-value-bytes) in
[`schemas.py`](../src/structcooker/schemas.py), then hands the config to
`datacooker pipeline`, which:

1. **sizes** every item — source-file `st_size` (build) or stored value bytes from
   `<db>.index.tsv` (rebuild);
2. **distributes** items across node shards by input size; each Ray driver respects
   allocated CPUs and adjusts task admission using measured memory pressure;
3. **submits** the job array → **merges** the shards → writes the target's
   own `.index.tsv` + `.meta.json`, all `afterok`-chained.

OOM retries and result buffering are bounded. Admission leaves memory headroom,
but cannot guarantee that an arbitrarily large individual record fits. Each run
has unique shard paths and an explicit merge manifest. Dependent planning waits
for successful upstream completion before reading generated inputs.

No human tunes n_jobs / mem / chunk / shards. The `schema:` tag never reaches the
engine: `structcooker` materializes a schema-stripped `engine.yaml` in the workdir.

## Op is inferred

| config has            | op        | source sizing            |
|-----------------------|-----------|--------------------------|
| `env_path` + `data_dir`/`file_pattern` | `build`   | source-file `st_size`    |
| `old_env_path` + `new_env_path`        | `rebuild` | `<source>.index.tsv`     |

## Schemas (see `schemas.py`)

| schema | value                              | key           | E     |
|--------|------------------------------------|---------------|-------|
| A      | `{assembly_dict, metadata_dict}`   | pdbid         | 100   |
| B      | cifmol_attached (per assembly)     | pdbid         | 100   |
| C      | chain BioMol                       | pdbid_chain   | 100   |
| D      | `{template_mols}`                  | pdbid_chain   | 20    |
| E      | `{msa_dict}`                       | seqid         | 100   |
| F      | raw scalar bytes                   | mixed         | 1     |
| G      | `{chem_comp_dict}`                 | comp_id       | 5     |
| H      | pruned CIFMol                      | pdbid_assembly_model_altloc | 100 |
| I      | `{assembly: {cifmol_dict: ...}}`   | pdbid         | 100   |

## Targets

DB destinations default to **`/data/shared/cssb_data/BioMol_clean/`** and can be
changed with `OUTPUT_ROOT`. Production `BioMol/` stays read-only throughout the
recovery, including after verification. Historical CCD and chain-selection
references under the output tree are provided reproduction inputs; preserve them.

## Layout

```
db/
  pdb/     cif.yaml  cif_attached.yaml          # A ingest, B attach
  train/   train_20210930.yaml  train_20260301.yaml   # B release-date filters
  chain/   …   template/   …   msa/   …   ccd/   …
  distillation/   …        # OpenFold distillation (OFD): short / long / rna / disordered
  teddymer/       cif  cif_attached  msa_search  msa  msa_d2k  msa_wo_lower  hmmsearch
  afdb_multimer/  {homo,hetero}dimer_cif(_attached)  msa  msa_d2k
  metadata/       …        # fasta / seq_id_map (shared) / <db>_seq_cluster30 (per DB)
```
