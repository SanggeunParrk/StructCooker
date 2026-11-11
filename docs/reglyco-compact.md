# Re-Glyco, stored compactly

The RFD4-Proteina Re-Glyco set (`train_datasets/streaming-latest/synthetic/glycans/2026_09_27_reglyco`, see
[raw materials](raw-materials.md)) is GlycoShape glycans rebuilt on PDB glycosylation sites: **2,992,448 records** of
11,160 PDB entries, **1,573.6 GB** as shipped (10.85 TB decompressed). This page records why that is mostly repeated
protein, how small a random-access copy can be, and the layout to build. Everything below is measured on one shipped shard
(`rs_00005.bin`, 700 records of 5 entries) and the Hub's two tables, `blob/index.parquet` and `examples.parquet`; the code
is `scripts/analysis/reglyco_compact.py`, `reglyco_extrapolate.py` (+ `.sbatch`), outputs in
`$BIOMOL_ROOT/materials/intermediate/reglyco_compact/`.

**Result.** One zstd delta per record against a base record of the same entry stores the whole set in **about 8 GB**
(7-9 GB measured range; 6.5 GB for the training view alone), 190 times smaller than shipped, with a record decoded in
0.1 ms and byte-identical to the original. Converting the same records to our LMDB format would not do this: the
records are already zstd-compressed one by one, and re-compressing them one by one gives 164 MB for the shard's 159 MB.

## What differs between records

A record is one annotated mmCIF named `<xx>/<pdb>/<pdb>_<site>_<glycan>_conf<NNN>.cif`: one PDB entry, one glycosylation
site, one GlycoShape glycan, one conformer. The index gives, for every record, `shard`, `offset`, `length`,
`uncompressed_length` and the sha256 of the decompressed CIF.

- The whole protein context is rewritten into every record: 14,376 lines for `2pmv_A293_GS00178_conf000`.
- Two conformers of one site differ in **83 of those lines**, all glycan `HETATM` rows (coordinates only).
- Another site or glycan of the same entry replaces the glycan rows: a group's first record costs 2.7 KB against the entry
  base, a further conformer 1.9 KB against its group's first.
- The population: 11,160 entries, 499,626 (entry, site, glycan) groups, 6.0 conformers per group, 268 records per entry
  (up to 3,261). 11,064 entries lie entirely in one shard. Records are 3.8 MB at the median, 9.4 MB at the 99th percentile
  and 13.7 MB at most.

## Layouts measured

700 records, 5 entries, 108 groups. Shipped: 159.23 MB (227 KB per record); decompressed 972.6 MB. Every layout below was
decoded again and compared with the original bytes and the index's `cif_sha256`: all 700 records match (A3: 31 records).

| Layout | What is stored | Size | Per record (A1/A2: the delta part) | Decode | Random access |
|---|---|---:|---:|---|---|
| shipped | one zstd frame per record | 159.2 MB | 227 KB | - | record |
| per-record, level 6 | the same, recompressed (what our `to_bytes` does) | 164.1 MB | 234 KB | - | record |
| B, level 12 | one solid zstd frame per group | 26.1 MB | 37 KB | 2 ms per group (9 MB) | group |
| A1, level 9 | entry base + one delta per record against it | 2.84 MB | 2.6 KB | 0.1 ms (+2 ms for a new entry) | record |
| **A2, level 9** | entry base + group base (delta vs entry base) + conformer deltas vs the group base | **2.31 MB** | 1.7 KB per conformer | 0.1 ms | record |
| A2, level 19 | as A2 | 1.99 MB | 1.3 KB per conformer | 0.1 ms | record |
| A3 (sample) | entry base + explicit line diff (copy / literal ops), zstd 19 | - | 2.1 KB | 1 ms | record |

- A "delta" is `zstd` with the base record as a raw-content dictionary (`ZstdCompressionDict(..., DICT_TYPE_RAWCONTENT)`,
  window sized to cover base + record). The decoder needs the base bytes and nothing else, and any zstd binding can load a
  raw-content dictionary.
- **Levels 6, 9 and 12 give the same size; level 3 does not work** (51.9 MB: the match finder does not find the 1.9 MB
  base). The delta layouts must be built at level 6 or higher; level 19 saves another 12-17 % for about ten times the encode
  time (A2: 8.6 ms per record at level 9, about 100 ms at level 19).
- B is 11 times worse than A2: it removes the repetition inside a group, not the protein shared by all groups of an entry.
- A3, a hand-made line diff, lands where A1 at level 19 lands (2.1 KB per record). The 80 changed glycan rows of coordinates
  are the information floor of the text form, about 2 KB; zstd with a dictionary reaches it without any CIF-specific code.
- Not measured: a structure-of-arrays copy (coordinates of the changed atoms as float32, about 1 KB per record). It would
  halve the deltas again but loses the byte-exact CIF and ties the store to one parser. Worth it only if the loader reads
  arrays and never the CIF text.

## Recommendation: A2, zstd level 9

Why A2 over A1: group bases cost 2.7 KB, conformers 1.9 KB, and a conformer's neighbour in the same group is the closest
base there is. Why not B: 11 times larger and no per-record access. Decode is not the cost to worry about; parsing the 1-4 MB
CIF the loader gets back is orders of magnitude slower than the 0.1 ms delta decode.

**Records.** Every Hub record keeps its key, the file name without `.cif` (`2pmv_A293_GS00178_conf000`). The base of a group
is its lowest conformer, the base of an entry is the first record by path. A value is
`kind u8 | n u32 | ref_len u8 | ref | payload`: kind 0 = plain zstd 19 (the entry base, `ref` empty), kind 1 = delta
(`ref` = the key of its base: the entry base for a group base, the group base for a conformer). Decoding is at most two
steps, and a loader keeps the entry base (and the current group base) of every entry it is reading in a small cache.

## Size for the whole set

Exact population from the two Hub tables; sizes extrapolated from the sample (`reglyco_extrapolate.py`). An entry's base is
modelled from the real shipped length of its first record (the measured ratio is 0.83, range 0.82-0.84); the delta means
come from the sample, the range from its smallest and largest per-entry mean.

| | records | entries | groups | shipped | A2, level 9 | A2, level 19 | A1, level 9 |
|---|---:|---:|---:|---:|---:|---:|---:|
| all | 2,992,448 | 11,160 | 499,626 | 1,573.6 GB | **8.4 GB** (7.2-9.4) | 7.0 GB (6.2-7.8) | 10.6 GB (9.7-12.3) |
| training view | 2,307,945 | 9,215 | 384,598 | 1,192.1 GB | **6.5 GB** (5.6-7.3) | 5.5 GB (4.8-6.0) | 8.3 GB (7.5-9.5) |

The training view is the one the configs use: deposited before 2024-01-01, `glycan_experiment_split == 'train'`, an
in-vocabulary glycan (2,307,945 rows, as the README says). Of the 8.4 GB, 2.6 GB are entry bases, 1.3 GB group bases and
4.6 GB conformer deltas. Add 5-10 % for LMDB pages and keys.

Build cost, from the sample's timing: reading the shards (10.85 TB decompressed) about 3 core-hours, encoding at level 9
about 7 core-hours, so on the order of ten core-hours (level 19: about 85).

**What the sample cannot tell.** Five entries are 0.04 % of the set, and their records are smaller than the population's
(1.4 MB against 3.6 MB on average). The base part scales with the real sizes; the delta part is measured per glycan, not
per protein. If deltas grew with the record size (they should not: the changed lines are the glycan's), the upper bound is
about 18 GB. The range in the table is the per-entry spread, not a confidence interval.

## How it is built here

One stage next to the other glycan/PDB-derived DBs, to be built after the download passed `verify --sha`:

- `db_configs/glycan/reglyco.yaml` (op `build`, new schema **K**, codec raw: the values are compressed already) writes
  `lmdb/glycan/reglyco_pdb.lmdb`. `examples.parquet` is copied as `metadata/reglyco_examples.parquet`; the sampler joins on
  `path` (the loader filters `train` / in-vocabulary / date and weights by `sampling_weight` there, not in the DB).
- Work items are the 499,626 groups plus the 11,160 entry bases, from a `materialize` stage over `blob/index.parquet`
  (`reglyco_groups.tsv`: key, base key, and for each record `shard`, `offset`, `length`, `cif_sha256`). A reader
  (`readers/reglyco.py`) reads a group's frames through the index (never a shard sequentially: shards hold bytes the index
  does not reference), a recipe (`workflows/ingest/reglyco_compact.py`) encodes the deltas, decodes each once and checks the
  original `cif_sha256`, and the build fails on any mismatch. An item holds one group plus its entry base in memory (at most
  a few conformers of about 4 MB).
- `_v_K` in `build/schemas.py` checks `kind`, `n`, `ref` and the chain depth. The release check samples records and decodes
  them against the index hash.
- The raw shards stay where the download put them (`materials/raw/rfd4_proteina`, exactly as downloaded). Deleting them
  frees 1.57 TB, but the compact DB could then only be rebuilt by downloading the set again.

## Open points

- The delta size of large entries and of entries with many deposited glycans is unmeasured; run
  `LEVELS=9 sbatch scripts/analysis/reglyco_compact.sbatch <shard>` on a few shards of other size before fixing the schema.
- Whether to keep the 23 % of the records outside the training view (holdout, pre-2024 failures): they cost about 1.9 GB;
  keeping them leaves the holdout available.
- The loader's CIF parser (AtomWorks `annotated_cif`) is the real decode cost; parsing each entry's protein once and storing
  arrays would be a further step, outside this page.

## Reproduce

```bash
source scripts/cluster/cssb.env
sbatch scripts/analysis/reglyco_compact.sbatch                 # layouts on rs_00005 -> compact_rs_00005.json (3 min)
sbatch scripts/analysis/reglyco_extrapolate.sbatch             # population + sizes -> extrapolation.json (10 s)
```

`examples.parquet` was fetched separately to `materials/intermediate/reglyco_compact/hub/`; the extrapolation reads it
from there.
