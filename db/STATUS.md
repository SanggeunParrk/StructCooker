# db/ migration ledger

Historical migration notes below do not establish current output correctness.
For the September 10 recovery, production comparisons, and final build status,
see [the recovery report](../docs/recovery-2026-09-10.md).

Every non-distillation BioMol DB and the status of its declarative `db/*.yaml`.
CIF reproduction also requires the explicit historical CCD reference described in
[the recovery report](../docs/recovery-2026-09-10.md); a missing or incompatible
definition raises an error instead of dropping observed atoms.
`structcooker build <name>` works only for ✅ rows; 🔴 rows need engine work first
(a build whose items are **keys** can't be sized by input `st_size`, and script-only
builds must be absorbed into a datacooker recipe/op before they get a config).

## ✅ Clean — config written, planning-first works

| db config | op | schema | source of truth |
|---|---|---|---|
| `pdb/cif` | build | A | configured raw CIF tree (real files → st_size) |
| `pdb/cif_attached` | rebuild | B | cif_pdb `.index.tsv` |
| `train/train_20210930` | rebuild | B | cif_pdb_attached `.index.tsv` |
| `train/train_20260301` | rebuild | B | cif_pdb_attached `.index.tsv` |
| `train/train_20210930_mpnn` | rebuild | H | train_20210930 `.index.tsv` |
| `msa/a3m` | build | E | real `*.a3m` files → st_size |
| `msa/a3m_d16k` / `a3m_d2k` / `a3m_d512` | rebuild | E | a3m `.index.tsv` (depth cap 16000/2000/512) |
| `ccd/ccd` | build | G | component `.cif` files → st_size |
| `chain/cif_chain` | rebuild | C | cif_pdb `.index.tsv` (explode per chain; optional production model-choice reference) |

Replaces the 10 hand-tuned `cif_pdb_{light,medium,large,xlarge,huge,monster,smoke,…}`
variants with `pdb/cif.yaml` alone. Schema **G** (CCD chem-component) added for `ccd`.

Schema **H** (pruned CIFMol, one key per assembly-model-altloc) added for the first
*model-input* DB, `train/train_20210930_mpnn` — a view of a CIF DB carrying only the
features one model reads. It is the one rebuild that is 1:N rather than 1:1
(`explode_entries`), so a key is a training sample. See `docs/mpnn_pruning.md` for the
size measurements and for why its `index_table` must be read through
`instructions.transforms.pruning.load_pruned_cifmol`.

### MSA depth-cap bypass absorbed (was a 🔴 row)
`scripts/maintenance/lightweight_msa.py` + `build_msa_long_d2k.sh` are replaced by the
`workflows/exports/cap_msa_depth.py` rebuild recipe (+ `adapt_msa_for_cap` adapter):
`db/msa/a3m_{d16k,d2k,d512}.yaml` differ only in `parameters.max_depth`. Verified on
a3m_d16k (178,249 keys, key-set complete, deep records capped to 16000, standard-frame
output). The recipe caps `msa_dict` in place; the adapter renames the input so the
same-named output does not silently return the *uncapped* record.

### Ray build/rebuild execution

Build/rebuild fan-out uses Ray. The driver admits work gradually using measured
memory pressure and respects `n_jobs`, CPU affinity, and the SLURM allocation.
Metadata is put once and installed once per worker per run; Python dictionaries
are deserialized in each worker, not shared as zero-copy dictionaries. Compact
SignalP metadata avoids distributing the full sequence lookup.

OOM retries are bounded. This is not a guarantee that every record fits in memory.
Input sizes are distributed once across shards. Result buffering is bounded by
count, bytes, and time. Each pipeline run uses unique shard paths and an explicit
merge manifest; dependent planning waits for successful upstream completion.

Run Ray tests and real builds on compute nodes, never the login node. The frozen
September 10 release passed 112 tests (SLURM job 266534). Older MSA measurements
apply to their historical implementation, not this recovery's performance claim.

### Superseded old configs deleted (git-recoverable)
`configs/ingest/cif_pdb_{light,medium,large,large_fast,xlarge,xlarge_fast,huge,monster}`,
`cif_lmdb`, `cif_lmdb_refactor_test`, `filters/train_filter{,_20260224}`, `ingest/a3m_lmdb`,
`ingest/ccd_lmdb` (14). Untracked-and-superseded, left for now (no git safety net):
`cif_pdb_large_retry`, `cif_pdb_smoke`, `metadata/attach_cifpdb_revisit`, `ingest/ccd_lmdb_20260718`.
`logs/` cleared (5894 files / 393 MB).

## 🟡 Multi-stage — belongs to the step-4 pipeline DAG, not a single config

| db | why |
|---|---|
| `valid/valid1` (stage 1) → attach → `valid/valid2` (stage 2) | stage 2 dedups sequences vs **train/valid_1 fasta + seq clusters** — cross-set deps (`validation_stage2*`), needs the orchestration DAG, not a lone rebuild config |

## 🟢 Metadata + projection ops (new op class)

`structcooker build` now auto-infers **materialize** (recipe → TSV/fasta) and
**extract** (LMDB → TSV) ops alongside build/rebuild, routed to
`datacooker.cli.workflow` as a single afterok-chainable job. The PDB metadata chain
is ported under `db/metadata/`:

| config | op | notes |
|---|---|---|
| `cif_fasta` / `cif_metadata` | extract | cif_pdb → fasta / metadata TSV (recipes verified on prod records) |
| `seq_id_map` | materialize | **optional seed exception** — seq_id is an assigned counter, so from-scratch ≠ production; seed via the provided reference `DATA_ROOT/reference/seq_id_map.tsv` to match, else fresh ids. Validated. |
| `{pdb,ofd,teddymer,afm}_seq_cluster30` | materialize | mmseqs2 at 30% (antibodies via cd-hit), **per DB** over that DB's own fasta; ids `c<DB>_<rep seq_id>` ([scheme](../docs/seq-id-and-cluster-scheme.md)). Replaced the global `seq_cluster40` / `seq_cluster30` (2026-09-23). SabDab is an external input |
| `interacting_seq_ids` / `interacting_seq_clusters` | extract / materialize | interface partners for valid-2 dedup |

### Key-list sizing gap — CLOSED
A build with `keyed: true` sizes items by **count** (uniform), not input `st_size`,
so key-list builds (items are seqids/chain-ids, resolved by lookup) plan cleanly.

### Template hmm pipeline — ABSORBED (4th op: `parallel`)
The full template pipeline is now on the clean surface, wired in MANIFEST order
`hmmsearch → template_candidates → template_seqs → seqid_template_mols → chain_template`:

| config | op | notes |
|---|---|---|
| `template/hmmsearch` | parallel | Phase 0 — hmmbuild(query a3m→HMM) + hmmsearch(vs L-chain DB). `split_recipe` → `datacooker.cli.workflow parallel-run`. Idempotent |
| `template/template_candidates` | materialize | Phase 1/2 — seqid_to_chains, chain_to_templates, reduced_hmm (verbatim port; logic verified on synthetic inputs) |
| `template/template_seqs` | materialize | Phase 3 seq maps (seqid_to_seq, chain_to_seq) — supersedes the old cif_chain_seq LMDB |
| `template/seqid_template_mols` / `chain_template` | build (keyed) | Phase 3/4, schema D |

`msa/msa_rna`, `valid/valid1`, `valid/valid1_attach`, `valid/valid2` are also ported.

### Minor derivations — DONE
The small list/fasta projections feeding the hmm pipeline are ported as materialize/
parallel ops: `metadata/pdb_polypeptide_L` (L-chain fasta = hmmsearch template DB),
`metadata/template_chain_filelist` (protein-chain work list), `metadata/seqid_template_filelist`
(Phase 3 keyed-build item list), `template/msa_wo_lower` (a3m insertion-strip). The full
**33-node MANIFEST** topo-sorts end to end.

## PDB reproduction — configuration coverage

Every step from raw mmCIF to the training/validation DBs is a `db/**` config; nothing
runs off a legacy script. What is *not* reproduced here (by design, provided like the
raw mmCIF/CCD downloads): the external-tool / raw inputs — **SabDab** (antibody summary),
**SignalP** outputs, and the raw per-chain **a3m** MSAs. The legacy `rekey_seq_id_db.py`
is a one-off width migration (`P0000007` → `P…020d`) for *pre-existing* DBs; a
from-scratch build already emits 20-width ids (`_SEQ_ID_WIDTH`), so it is not a
reproduction step.

### Historical manual mmCIF fixes

The source session recorded 53 PDB entries (`db/pdb/manual_cif_fixes.txt`) that failed
with its wwPDB snapshot — mostly NMR ensembles whose non-polymer ligand is re-numbered per model,
breaking the atom→scheme match (verified on `1ai0`: `IPH` at `auth_seq` 22 in some models,
31 in others, vs a single scheme number). The production build substituted an older
known-good cif for each before ingest (legacy `scripts/manually_fix_cif.py`, source
`BioMolDB_2024Oct21`). Ported as `structcooker fix-cif --source <corrected-cif-dir>`; run
it for a snapshot requiring those substitutions, and `inspect` reports the state.
These historical substitutions are distinct from the CCD compatibility corrections
verified in the September 10 recovery.
Full write-up: [docs/manual-cif-fixes.md](../docs/manual-cif-fixes.md).
**Not needed with the current parser (checked 2026-10-01):** all 53 records were compared with
production `BioMol/lmdb/pdb/cif/cif_pdb.lmdb` -- the 48 that build are identical in every array,
and the other 5 (2g10, 2icy, 2q44, 4xq2, 9gdy) are absent from production too.

**Superseded run:** the 233,579-entry database described above used the old input
snapshot. The current recovery reads `BioMol/materials/raw/cif` and has 249,676 CIF
keys, with no missing production keys. It does not wait for the old manual-fix
snapshot. See the recovery report for the actual comparison scope, CCD compatibility
references, deliberate source differences, and final build status.

The OpenFold3 distillation sets. **Recipes exist** (`workflows/ingest/openfold_*`);
configs are ported onto the clean `db/distillation/` surface (env-var paths, schema tag,
ops knobs dropped — same shape as long_cif). Production stores the *downstream* DBs
(`cif_*_attached`, capped `msa_*_{d2k,d512,d5k}`, `template_topn_*`/`_disordered`), so the
full chain is base → derived.

**✅ ported + resolve + small-validated (21 configs):**

| kind | schema | configs |
|---|---|---|
| structure (base cif) | A | `{long,short,rna,disordered}_cif` |
| msa (base, depth 65536) | E | `{long,short,rna,disordered}_msa` |
| msa depth-caps | E | `{long,disordered}_msa_{d2k,d512}`, `{short,rna}_msa_{d2k,d512,d5k}` |
| template (top-N) | D | `{long,short}_template_topn` |
| template (disordered) | D | `disordered_template` |

Base-structure pipeline validated in-process (Ray-free): one real `structure.npz`
(`MGYP003648360693`) → recipe → schema-A valid → serialize/deserialize round-trips.

**✅ `cif_*_attached` (schema B) — ported (25 configs total).** Previously blocked on a
rewrap step + missing metadata; both resolved. The base `cif_{set}` already emits the
wrapped `assembly_dict`/`metadata_dict` layout (the openfold_structure recipe), so no
separate rewrap is needed, and `seq_id_map` / the per-DB clusterings are now on the clean
surface (see the metadata section). `cif_{long,short,rna,disordered}_attached` are plain
A→B rebuilds of the base + the shared attach recipe (proven by pdb/cif_attached).

The earlier sample-only restriction was superseded by the explicit end-to-end
distillation build request. See [the September 15 follow-up](../docs/distillation-followup-2026-09-15.md)
for the completed long CIF audit and active native SLURM MSA/template graph.
Submission is not completion; terminal coverage and validation reports decide status.

## Release status — 2026-10-07

End to end: `structcooker build-all --manifest db/MANIFEST_all.yaml` (94 stages; downloads,
staging and reference inputs in [docs/e2e-build.md](../docs/e2e-build.md)).

| set | cif / cif_attached (structures) | msa (+ caps) (unique sequences) | template (unique sequences with ≥ 1 candidate) | b_factor |
|---|---|---|---|---|
| PDB | ✅ 249,676 / 249,644 (`cPDB_`) | ✅ 178,249 (+d16k/d2k/d512); RNA 6,572 | production DB (inputs unchanged) | experimental |
| OFD short | ✅ 430,418 / 430,418 (released relaxed PDB) | ✅ 430,245 (+d2k/d512/d5k) | ✅ 430,370 (rebuilt on the fixed cif_chain) | pLDDT |
| OFD long | ✅ 16,099,404 / 16,099,404 (released relaxed PDB; 80 ids have no model in the release) | ✅ 16,098,794 (+d2k/d512) | ✅ 16,080,773 (1,644,360 with hits on the added chains patched) | pLDDT |
| OFD rna | ✅ 126,778 / 126,778 | ✅ 126,751 (+d2k/d512/d5k) | — | 0 (not published) |
| OFD disordered | ✅ 28,567 / 28,567 (+contacts) | ✅ 19,649 (+d2k/d512) | ✅ 78,176 | 0 (not published) |
| teddymer | ✅ 510,454 / 510,454 | ✅ 999,853 (+d2k) | ✅ 998,878 | pLDDT |
| AFDB homodimer | ✅ 1,750,755 / 1,750,755 | ✅ 1,721,635 (+d2k, shared) | ✅ 1,720,477 (shared, by seq_id) | pLDDT |
| AFDB heterodimer | ✅ 80,248 / 80,248 | (above) | (above) | pLDDT |
| train / valid | ✅ train_20210930 167,912, train_20260301 233,579, mpnn 484,383; valid1 21,527, valid2 617 | | | |

`chain/cif_chain` gained 293 PDB entries (51,895 chains) on 2026-09-30: an atom cap dropped a
whole entry when any assembly was over it. Template DBs built before that miss those chains as
candidates; the OFD short and long templates were rebuilt/patched for that on 2026-10-06/07.
Partial rebuilds use `structcooker patch NAME --keys|--files [--nodes N]`.

MSA and template DBs are keyed by seq_id, so they count unique sequences, not structures: the
AFM release predicts one model per UniProt entry, and 126,492 homodimers share their sequence
with another entry (other strains, or proteins conserved across close species), so 1,750,755 +
80,248 structures hold 1,721,635 sequences. Template DBs hold the sequences with at least one
candidate released by 2021-09-30 (AFM: 1,159 have none).

## Teddymer and AFDB multimer — built 2026-09-28 (history)

| set | cif | cif_attached | msa / msa_d2k | template |
|---|---|---|---|---|
| teddymer (`MANIFEST_teddymer`) | ✅ 510,454 | ✅ 510,454 | ✅ 999,853 / 999,853 | ✅ 998,878 |
| AFDB homodimer (`MANIFEST_afdb_multimer`) | ✅ 1,750,755 | ✅ 1,750,755 | ✅ 1,721,635 / 1,721,635 (shared with heterodimer, keyed by seq_id) | ✅ 1,720,477 |
| AFDB heterodimer | ✅ 80,248 | ✅ 80,248 | (above) | (above) |

Raw inputs and their selection: [docs/raw-materials.md](../docs/raw-materials.md).

