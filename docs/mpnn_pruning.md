# Pruning a CIF DB into a model-input DB

`db/train/train_20210930_mpnn.yaml` derives an MPNN training DB from
`cif_pdb_attached_train_20210930_res9_chain300`. It is the first *model-input*
DB — a view of the general CIF DB carrying only what one model reads — so this
note records what the measurements were, in case the next one wants the same
treatment.

Everything below was measured by re-serializing real records from the source DB
at the same zstd level, not estimated.

## What the source costs

167,912 keys / 47.9 GiB. Each key holds **2.90 sub-entries** on average, one per
`{assembly}_{model}_{altloc}`. That multiplicity is mostly altlocs (46.8% of
records have >1 altloc, usually 3), not assemblies (31.7% have >1); extra models
are effectively absent (0.1%).

Per record, uncompressed, the largest fields:

| field | KB/rec | note |
|---|---:|---|
| `atoms.nodes.model_xyz` | 1916 | CCD *reference conformer* coords, kept as `<U7` strings — not NMR models |
| `atoms.edges.bond_*` | 1587 | CCD bond graph (value + 3× src/dst pairs) |
| `atoms.nodes.xyz` | 548 | float64 |
| `atoms.nodes.id` | 365 | `<U4` |
| `index_table.*` | 264 | **JSON integer lists in the header**, not arrays |
| everything else | 2320 | b-factor, occupancy, element, charge, residue labels, … |

Two of those are worth flagging because they are easy to read wrong:

* `model_xyz` is `chem_comp_atom.model_Cartn_x/y/z` attached from the CCD
  (`cif_instructions.py`), i.e. idealized template coordinates. It is not a
  multi-model NMR ensemble.
* `index_table` looks free if you only sum the array payloads, because
  `IndexTable.to_dict()` calls `.tolist()` and the arrays land in the JSON
  header instead. It is not free.

## What each pruning step buys

Cumulative, measured on 300 random records, as a fraction of the 47.9 GiB source:

| layout | vs source | DB |
|---|---:|---:|
| keep-list only (dtypes and index table unchanged) | 47.9% | 23.0 G |
| + `index_table` as int32 arrays, CSR dropped | 45.0% | 21.6 G |
| + `xyz` float64 → float32 | 38.0% | 18.2 G |
| + backbone N/CA/C/O atoms only | 17.1% | 8.2 G |

The headline: **field pruning alone stops at ~48%, not the ~15% a raw-byte
share suggests.** What gets dropped (CCD strings, bond graphs, `<U`-dtype
labels) is exactly what zstd was already compressing 20-40×, while float64
coordinates barely compress. Judge a pruning plan on compressed bytes.

## What this DB actually does

Only the first two rows: **keep-list + compact index table, full all-atom
float64 `xyz`**, plus the key explode. Backbone-only filtering and the float32
downcast are deliberately *not* applied — they change what the data means
(every nucleic-acid and ligand chain loses all its atoms, and 1.5% of records
end up with no atoms at all), and that is a modelling decision, not a storage
one.

### The explode costs disk

Splitting each record's sub-entries into their own keys removes the compression
window zstd was using to deduplicate near-identical altloc copies:

| | keys | DB |
|---|---:|---:|
| pruned, nested per pdb id | 167,912 | 21.6 G |
| pruned, exploded per sub-entry | ~517,000 | **~30 G** |

Paid for on purpose: one key is one training sample, so a read decompresses
~80 KB instead of ~250 KB, and the dataloader indexes/shuffles/shards on keys
rather than on (key, sub-entry) pairs.

If the extra 9 G ever matters, the lever is not compression — it is deciding
that near-identical altlocs should not each be a separate training sample.
Keeping one altloc per (assembly, model) halves the DB to ~18 G; assembly 1
only takes it to ~12.6 G.

## Reading it

The stored `index_table` has only `atom_to_res` and `res_to_chain`; the four CSR
arrays are a pure function of those two and are rebuilt on load. So
`CIFMolAttached.from_dict(value)` will **not** work directly — go through:

```python
from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.instructions.transforms.pruning import load_pruned_cifmol

mol = load_pruned_cifmol(from_bytes(raw_value))     # -> CIFMolAttached
```

`restore_index_table` is a no-op on records that already carry the full index
table, so a loader can call it unconditionally.

## Engine change this needed

`rebuild_lmdb` was 1:1 — one source key, one destination record. `explode_entries`
makes it 1:N over a `split_entries` result (`datacooker/lmdb/core.py`,
`tests/test_lmdb_explode.py`). Sharding, merge, index and resume all still work:
shards partition by source key so exploded keys stay disjoint, and `skip_existing`
recovers a destination record's source key from its prefix — which is why a
source key must not contain the separator.
