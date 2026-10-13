# Metadata

Metadata workflows assign stable sequence ids and attach per-sequence metadata
onto structure records. The id and cluster conventions are described in the
[seq-id and cluster scheme](../seq-id-and-cluster-scheme.md).

## Sequence id map

Assign one id to every unique sequence across **all** DBs, so a3m, templates,
attach, and clusters join on the same key.

- Recipe: `structcooker.workflows.metadata.seq_id_map` (`seq_id_map`)
- Config (materialize):
  [`db/metadata/seq_id_map.yaml`](https://github.com/SanggeunParrk/StructCooker/blob/main/db/metadata/seq_id_map.yaml)
- Writer: `...writers.projections.write_seq_id_map` → `OUTPUT_ROOT/metadata/seq_id_map.tsv`
- Inputs (`fasta_paths`, in order): `cif_pdb.fasta`, the four OpenFold FASTAs
  (short / long / rna / disordered), `teddymer.fasta`, `afm.fasta`,
  `ofd_long_recovered.fasta`

The seq id is an assigned counter (mol-type prefix + ordinal), not a hash, so
the order of `fasta_paths` matters: a DB appended at the end leaves every
earlier id unchanged. `old_seq_id_map_path` (`DATA_ROOT/reference/seq_id_map.tsv`)
seeds the map: present → ids match production; absent → a fresh, self-consistent
id space. `structcooker inspect` reports which case applies.

```bash
pixi run structcooker build metadata/seq_id_map
```

Per-DB views of the shared map are separate materialize configs:
`metadata/teddymer_seq_id_map`, `metadata/afm_seq_id_map`
(`metadata.subset_seq_id_map`), and `metadata/distillation_monomer_seqid`
(`metadata.monomer_seqid`, entry → seq id for the OpenFold MSA keys).

## Attach sequence metadata

Rebuild a base CIF LMDB (schema `A`) into an *attached* LMDB (schema `B`),
enriching each `CIFMol` with sequence-level metadata (seq ids, clusters, ...).
This is a `rebuild`: `old_env_path` is read through a `reader.adapter`, a
`metadata_recipe` loads shared metadata once, `split_entries: true` fans each
entry into its assemblies, and records are written to `new_env_path`.

- Recipe: `structcooker.workflows.metadata.attach` (`cifmol_attached_dict`)
- Metadata recipe: `metadata.load_sequence_metadata`
  (`seq_metadata_map`, `seqid2seq`, `seqclusters2seqids`)
- Reader adapter: `...transforms.cifmol.convert_to_cifmol_dict`
  (`convert_to_cifmol_transformed` for the already-filtered valid-1 DB)
- `metadata_input`: `raw_fasta_path`, `seqid2seq_path` (`seq_id_map.tsv`),
  `seqcluster_path` (that DB's `seq_cluster30_<DB>.tsv`), `db_code`

| Config | Source → target | `db_code` |
| --- | --- | --- |
| `db/pdb/cif_attached.yaml` | `lmdb/pdb/cif/cif_pdb.lmdb` → `lmdb/pdb/cif_attached/cif_pdb_attached.lmdb` | `PDB` |
| `db/distillation/{long,short,rna,disordered}_cif_attached.yaml` | `lmdb/distillation_<set>/cif/` → `.../cif_attached/` (disordered then gets contacts from `distillation/disordered_contacts`) | `OFD` |
| `db/teddymer/cif_attached.yaml` | teddymer CIF → attached | `TDM` |
| `db/afdb_multimer/{homo,hetero}dimer_cif_attached.yaml` | AFM CIF → attached | `AFM` |
| `db/valid/valid1_attach.yaml` | `cif_valid_1.lmdb` → `cif_attached_valid_1.lmdb` | `PDB` |

```bash
pixi run structcooker build pdb/cif_attached
```

## Extract metadata

`structcooker.workflows.metadata.extract` (`metadata_dict`) dumps per-chain
metadata (resolution, dates, chain/residue/atom counts, NA/D-form flags) from
`cif_pdb.lmdb` into `OUTPUT_ROOT/metadata/cif_metadata.tsv` for the template
date rule and downstream filtering.

- Config (extract): `db/metadata/cif_metadata.yaml`
  (`writer.materializer: ...projections.write_metadata`)

```bash
pixi run structcooker build metadata/cif_metadata
```
