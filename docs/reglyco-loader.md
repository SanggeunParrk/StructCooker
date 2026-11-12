# Reading the Re-Glyco DB (for loader authors)

`lmdb/glycan/reglyco/reglyco.lmdb` holds the RFD4-Proteina Re-Glyco set as **ordinary BioMol records**
(`{assembly_dict, metadata_dict}`, the shape of a `pdb/cif` value). It is the one DB of the release whose stored values are
not one plain `to_bytes` frame per key: the records of an entry are near-copies of each other, so all but the first are
stored as small deltas against a base record of the same entry ([why and how big](reglyco-compact.md)). A loader does not
need to know that beyond one rule: **read through `DeltaReader` / `read_record`, never `from_bytes` on a raw value.**

```python
from datacooker.lmdb.sharded import open_env
from structcooker.instructions.transforms.delta_codec import DeltaReader

with open_env("lmdb/glycan/reglyco/reglyco.lmdb", readonly=True, lock=False) as env:
    reader = DeltaReader(env)                            # keeps the base records it has decoded
    value = reader.read("2pmv_A293_GS00178_conf012")     # {"assembly_dict": ..., "metadata_dict": ...}
```

`value` is what `from_bytes` returns for any other schema-A DB: the same arrays, scalars and nesting, equal to the record
converted from the downloaded CIF. Everything after this line is the loader's business as for `cif_pdb`
(`convert_to_cifmol_dict`, cropping, featurisation).

## What to know

- **Keys** are the Hub's record names without `.cif`: `<pdb>_<site>_<glycan>_conf<NNN>` (`2pmv_A293_GS00178_conf012`). The
  keys of an entry share the `<pdb>_` prefix; `<site>_<glycan>` is a group (one site, one glycan; about 6 conformers).
- **Sampling metadata is not in the value.** The Hub's 90 per-record columns (glycan, site, UniProt, sampling weights, QC
  flags, `glycan_experiment_split`, `has_in_vocabulary_glycan`, size fields, ...) are in the sidecar
  `metadata/reglyco_examples.parquet`, one row per record, joined on its `key` column (the Hub's `path` column stays too;
  columns of the glycan view that differ from the full-assembly ones are prefixed `gv_`). `deposition_date` and
  `resolution` are also in `metadata_dict`. Filters of the RFD4 configs (deposition before 2024-01-01, split `train`,
  in-vocabulary glycan: 2,307,945 of 2,992,448 records) are applied on the sidecar, then the keys are read here.
- **Read the records of one entry together.** `DeltaReader` keeps the entry base and the group bases it has decoded
  (`cache_bytes`, 512 MB by default); the first read of an entry decodes its bases, the next ones reuse them. A loader that
  shuffles across all keys pays the base decode on most reads; sampling an entry's records in a batch, or keeping one
  `DeltaReader` per worker with a larger cache, avoids that. Measured on the built DB: a record costs about 51 ms to read with a cold reader and 36 ms with the entry bases cached (the
  time is mostly rebuilding the arrays); [the measurements](reglyco-compact.md#result-of-the-build).
- **A raw `txn.get` value is not a BioMol frame.** `codecs.from_bytes` on it raises a `ValueError` that names
  `read_record`. The entry bases (the first record of each entry by key) are plain frames and `from_bytes` reads them, but
  do not rely on which record is a base.
- **`read_record(env, key)`** is the one-shot form (no cache). `env` is whatever `open_env` returns (an LMDB environment or
  a shard collection).
- **What a record does not hold.** The AtomWorks per-atom columns (PLIP `condition_*` / `mask_*`, `nhyd`,
  `atom_array_index`), the imputation statistics, the CATH labels and the `_link_restore` note are not in the BioMol value;
  see [what the record does not hold](reglyco-compact.md#as-biomol-records). The CIFs they came from are in the raw
  download.
- **Patching and rebuilding** work per entry (a changed base invalidates its deltas); `structcooker patch` takes an entry's
  keys together.
