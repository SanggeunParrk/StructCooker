# BioMol release

Structure, sequence, MSA and template databases for protein / nucleic-acid structure
models, built from the PDB, the OpenFold3 distillation sets (OFD), Teddymer (TDM) and the
AlphaFold DB multimer set (AFM) by [StructCooker](https://github.com/SanggeunParrk/StructCooker).
Everything — downloads, supplied inputs, build intermediates and the databases — lives in
this one tree; the build resolves every path against it (`BIOMOL_ROOT`).

| page | |
|---|---|
| [index.html](index.html) | the release on one page: totals, data sets, what to copy, a searchable list of every database |
| [databases.md](databases.md) | every database: path, record schema, key, count, size, what it is built from |
| [scenarios.md](scenarios.md) | what to copy for a given use (PDB only, with MSA/templates, distillation, inference, rebuild, ...) |
| [provenance.md](provenance.md) | sources, versions, build dates, known gaps and differences |
| [check.md](check.md) | the latest `structcooker release-check` of this tree |
| [licenses.md](licenses.md) | the terms each source's data comes under |

## Layout

```
<root>/
├── docs/                   these pages (generated: structcooker release-docs)
├── lmdb/                   the databases, one LMDB per path (see databases.md)
│   ├── ccd/                chemical components: biomol_CCD (build), ccd_inference (inference)
│   ├── pdb/                cif, cif_attached (+ train sets), msa, valid
│   ├── chain/              per-chain structures (cif_chain) and the PDB template DBs
│   ├── distillation_{short,long,rna,disordered}/   OpenFold3 distillation sets
│   ├── teddymer/  afdb_multimer/
│   └── foldbench/          what FoldBench needs beyond the PDB-wide DBs
├── metadata/               released tables: seq_id_map, per-DB clusters, cif_metadata, ...
├── materials/
│   ├── raw/                downloads exactly as fetched (mmCIF, CCD, OpenFold, Teddymer, AFDB, SabDab)
│   ├── reference/          supplied inputs the downloads cannot reproduce (seq_id seed, PDB MSAs, ...)
│   └── intermediate/       anything derived on the way (fastas, search outputs, work lists, ...)
└── release_check.json
```

A large shard collection (`<db>.lmdb/shards.json`) keeps its shards inside its own
directory (`<db>.lmdb/shards/NNN.lmdb`); copy the whole `<db>.lmdb` directory. Every DB has
two sidecars next to it: `<db>.lmdb.index.tsv` (every key and its value size) and
`<db>.lmdb.meta.json` (schema, entry count, bytes).

## Reading a record

```python
from datacooker.lmdb.sharded import open_env                    # plain LMDB or shard collection
from structcooker.instructions.transforms.codecs import from_bytes

with open_env("lmdb/pdb/msa/a3m_d2k.lmdb", readonly=True, lock=False) as env, env.begin() as txn:
    record = from_bytes(txn.get(b"P00000000000000000373"))      # {"msa_dict": {...}}
```

Values are BioMol records compressed with zstd (`biomol_zstd`); `from_bytes` decodes them.
Without StructCooker installed, `lmdb` + `zstandard` + the BioMol codec are needed.

## Ids

- **seq_id** — one id space shared by every database: a type letter and 20 digits,
  `P00000000000000000373`. Letters: `P` protein (L), `Q` protein (D), `R` RNA, `D` DNA,
  `N` DNA/RNA hybrid, `B` branched (glycans), `L` ligand, `X` unknown. The same sequence
  has the same seq_id in every database, so cross-database joins go through it
  (`metadata/seq_id_map.tsv`: seq_id, sequence).
- **cluster id** — clustering is per database (30% identity, 80% coverage, MMseqs2):
  `c{DB}_{seq_id}` with DB = `PDB`, `OFD`, `TDM`, `AFM` (`metadata/seq_cluster30_<DB>.tsv`:
  cluster, members).
- **structure keys** — PDB id (`1abc`), assembly (`1abc_1_1_.`: id, assembly, model,
  altloc), chain (`1abc_A`), OpenFold entry (`MGYP…`, RNA accession), Teddymer
  (`<index>DI_AF-…_TED01_TED02`), AFDB entity (`AF-0000000065760022`).
- **MSA and template DBs** are keyed by seq_id (one record per distinct sequence);
  the PDB chain template DB by chain.

## Record schemas

| | value | key |
|---|---|---|
| A | structure: `assembly_dict` + `metadata_dict` | entry |
| B | structure with sequence metadata attached (seq_id, cluster) | entry |
| C | one chain | chain |
| D | templates (`template_mols`) | seq_id or chain |
| E | MSA (`msa_dict`: aligned sequences, deletions, profile, headers) | seq_id |
| G | CCD component (`chem_comp_dict`) | component code |
| H | pruned structure for ProteinMPNN | assembly/model/altloc |
| I | filtered structure (validation) | entry |
