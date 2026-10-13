# FoldBench inference DBs

FoldBench (BEAM-Labs, 9 tasks, 1,522 assemblies) is predicted from lookups by key. Every sequence has already been
searched and every template resolved for the whole PDB, so a FoldBench DB is a **key-list subset** of the full-release
LMDBs, not a new search. `structcooker foldbench` derives the key lists from a FoldBench package, copies exactly those
records byte-for-byte and verifies that every lookup a `data.yaml` makes resolves.

| asset | out path | key | value | cut from |
|---|---|---|---|---|
| protein MSA | `msa/a3m.lmdb` | `P` + 20 digits | `{msa_dict}` (schema E) | `lmdb/pdb/msa/a3m.lmdb` |
| RNA MSA | `msa/a3m_rna.lmdb` | `R` + 20 digits | `{msa_dict}` (schema E) | `lmdb/pdb/msa/msa_rna.lmdb` |
| template | `template/template.lmdb` | `{pdbid}_{chain}` | `{template_mols}` (schema D) | `lmdb/chain/template/template.lmdb` |
| CCD | `CCD/<name>.lmdb` | component code | CCDMol graph (`atoms` / `residues` / `chains`) | `lmdb/ccd/ccd_inference.lmdb` |

Sources are the release tree (`BIOMOL_ROOT`, default `/data/shared/cssb_data/BioMol_clean`), at the paths its db
configs write (`msa/a3m`, `msa/msa_rna`, `template/chain_template`). What FoldBench needs beyond those also lives in the
release's `lmdb/`, and `build` uses it by default:

- `lmdb/ccd/ccd_inference.lmdb` -- the CCD as inference reads it: one CCDMol graph per component, a pure reshape of
  the release's schema-G `lmdb/ccd/biomol_CCD.lmdb` (`scripts/staging/build_inference_ccd.py`; UNL, which has no atoms, is the
  one component left out). For FoldBench's 956 components it equals the July `CCD_20260718.lmdb` value for value.
- `lmdb/foldbench/msa_edge.lmdb` -- the four from-scratch edge MSAs (below).

## Package input

A FoldBench package is the directory `inference` reads from:

```
<package>/inputs/index.txt  inputs/<target>/data.yaml                # a3m_key, template per chain: the lookups to resolve
<package>/meta/chain_seqid_map.tsv                                   # target, chain_id, entity_type, seq_id, ..., sequence
<package>/meta/seqids_protein.txt  seqids_rna.txt                    # seq_ids; only used when inputs/ is absent
```

MSA and template keys are the `a3m_key` and `template` values of every `data.yaml` (so they include the from-scratch edge
ids); CCD codes are the inline `(XXX)` residues and ligands of `chain_seqid_map.tsv`.

## Build

```bash
PKG=/data/psk6950/foldbench_src                                                   # inputs/, meta/
structcooker foldbench keys   $PKG                                                 # counts
structcooker foldbench build  $PKG OUT                                            # SLURM job; --local to run here
structcooker foldbench verify $PKG OUT
```

`build` never replaces an existing DB, reports keys no source holds, and writes `OUT/foldbench_build.json` (requested,
written, from_fallback, missing, bytes per asset).

## Compared with the earlier FoldBench DBs

The earlier DBs (`/data/shared/cssb_data/BioMol/foldbench`) were cut from the July `BioMol` tree. Cut from the release
tree instead (checked 2026-10-13 on their own key lists): templates are byte-identical (1,937 found, 987 absent in
both); protein MSAs carry the same alignments, 32 of a 300-key sample differing only in the `headers/database` labels;
RNA MSAs differ in depth for most keys, as the release rebuilt them from the sanitized RNA alignments
([RNA MSA exceptions](history/rna_msa_sanitized_exceptions.md)).

## What is not rebuilt here

- **Edge sequences** (3 chromophore-abbreviated entries and a 99-mer, ids `P0300001`-`P0300004`) are not in the PDB-wide DB.
  They were searched from scratch once with `instructions.transforms.msa.run_msa_search` (SignalP, hhblits UniRef30 and BFD);
  They are kept in `lmdb/foldbench/msa_edge.lmdb`, which `build` uses as its fallback.
- **DNA** has no MSA (matching AF3): `D` / `N` / `Q` keys are expected to be absent and are counted separately by `verify`.
- **Templates** are limited by the cutoff the full-release template DB was built with (2021-09-30); this command does not
  re-filter by date.
