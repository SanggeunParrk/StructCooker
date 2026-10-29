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
PKG=$BIOMOL_ROOT/materials/reference/foldbench                                    # inputs/, meta/
structcooker foldbench keys   $PKG                                                 # counts
structcooker foldbench build  $PKG OUT                                            # SLURM job; --local to run here
structcooker foldbench verify $PKG OUT
```

The built DB lives in the release itself: `OUT = $BIOMOL_ROOT/lmdb/foldbench`, next to `msa_edge.lmdb`, so that
`lmdb/foldbench/{msa,template,CCD}` is the FoldBench DB on every server that holds the release tree (H100 and cssb:
`BioMolDB_20261014/lmdb/foldbench/`). The 2026-10-29 build there has 1,568 protein / 58 RNA MSAs (RNA tokens A U G C,
the release order), 1,941 of 2,395 requested templates (454 are absent from the release, all DNA/RNA chains: protein templates are only built for polypeptide chains) and 956 CCD components.

`build` never replaces an existing DB, reports keys no source holds, and writes `OUT/foldbench_build.json` (requested,
written, from_fallback, missing, bytes per asset).

## Compared with the earlier FoldBench DBs

The earlier DBs were cut from the July `BioMol` tree (since deleted). Cut from the release
tree instead (checked 2026-10-13 on their own key lists): templates are byte-identical (1,937 found, 987 absent in
both); protein MSAs carry the same alignments, 32 of a 300-key sample differing only in the `headers/database` labels;
RNA MSAs differ in depth for most keys, as the release rebuilt them from the sanitized RNA alignments
([RNA MSA exceptions](history/rna_msa_sanitized_exceptions.md)).

## What is not rebuilt here

- **Edge sequences** (3 chromophore-abbreviated entries and a 99-mer, ids `P0300001`-`P0300004`) are not in the PDB-wide DB.
  They were searched from scratch once (SignalP, then HHblits against UniRef30 and BFD; that search code has since been
  removed) and are kept in `lmdb/foldbench/msa_edge.lmdb`, which `build` uses as its fallback.
- **DNA** has no MSA (matching AF3): `D` / `N` / `Q` keys are expected to be absent and are counted separately by `verify`.
- **Templates** are limited by the cutoff the full-release template DB was built with (2021-09-30); this command does not
  re-filter by date. A found key can hold **zero** templates: the recipe (`db/template/template_candidates.yaml`) takes the top 20
  hmmsearch hits released at least 60 days before the query, then drops any hit whose alignment covers more than 95% of the
  query (`max_query_coverage`), so a chain whose every hit is a near-copy ends up with an empty `template_mols`. About 20% of
  of the original chains and 33% of the 17,356 added on 2026-10-29 are in that state.
- **Template history.** The first 2026-07 template build used a `seqid_to_chains.tsv` older than the chain list, so 17,356
  polypeptide chains (7,065 sequences, incl. 31 FoldBench chains) had no key. They were added on 2026-10-29 with
  `scripts/template_fix/` (phases 1-4 on that subset, then `merge_into_release.py`); the release `template.lmdb` is now
  1,025,284 chains and `seqid_template_mols.lmdb` 178,611 sequences. 118 sequences (1,614 chains, 97 with an MSA) never got an
  hmmsearch output in either build and still have no templates.
