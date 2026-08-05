# `structcooker build-all` — the BioMol DAG

One command turns [`db/MANIFEST.yaml`](../db/MANIFEST.yaml) into a datacooker
meta-recipe: **31 database configs**, topologically sorted, each dispatched by the op
class `structcooker build` infers from its config, and SLURM-chained with `afterok` so a
node only starts once its upstreams exit 0. Already-built outputs are skipped.

> A styled, interactive version of these diagrams is published as an artifact:
> <https://claude.ai/code/artifact/99280998-fe04-4500-be9f-9746ad8ef299>
> (re-render from [`docs/build-all-dag.html`](build-all-dag.html) if it needs updating).

- **31** db nodes · **5** op classes · **3** deliverables · `afterok` chaining
- Source of truth: `db/MANIFEST.yaml` · op inference: `src/structcooker/cli.py:_op_of`
- Run `structcooker inspect` before `structcooker build-all`.

## 1. The dependency graph

Left → right is build order. Each node is colored by its op. Dashed grey nodes are
external raw inputs — **not** DAG nodes; they are fetched by `structcooker download` /
provided, and checked by `structcooker inspect`.

```mermaid
graph LR
  ext_mmcif["mmCIF (wwPDB)"]:::external
  ext_a3m["a3m MSAs"]:::external
  ext_rna["RNA a3m"]:::external
  ext_sab["SabDab"]:::external
  ext_sig["SignalP out"]:::external

  ccd["ccd/ccd"]:::build
  a3m["msa/a3m"]:::build
  rna["msa/msa_rna"]:::build
  cif["pdb/cif"]:::build

  fasta["metadata/cif_fasta"]:::extract
  meta["metadata/cif_metadata"]:::extract
  seqid["metadata/seq_id_map"]:::materialize
  sc40["metadata/seq_cluster40"]:::materialize
  sc30["metadata/seq_cluster30"]:::materialize
  isi["metadata/interacting_seq_ids"]:::extract
  isc["metadata/interacting_seq_clusters"]:::materialize
  ppl["metadata/pdb_polypeptide_L"]:::materialize
  tcf["metadata/template_chain_filelist"]:::materialize
  stf["metadata/seqid_template_filelist"]:::materialize

  att["pdb/cif_attached"]:::rebuild
  chain["chain/cif_chain"]:::rebuild

  d16k["msa/a3m_d16k"]:::rebuild
  d2k["msa/a3m_d2k"]:::rebuild
  d512["msa/a3m_d512"]:::rebuild

  wol["template/msa_wo_lower"]:::parallel
  hmm["template/hmmsearch"]:::parallel
  tcand["template/template_candidates"]:::materialize
  tseq["template/template_seqs"]:::materialize
  stm["template/seqid_template_mols"]:::build
  ctpl["template/chain_template"]:::build

  tr09["train/train_20210930"]:::rebuild
  tr26["train/train_20260301"]:::rebuild
  mpnn["train/train_20210930_mpnn"]:::rebuild
  v1["valid/valid1"]:::rebuild
  v1a["valid/valid1_attach"]:::rebuild
  v2["valid/valid2"]:::rebuild

  ext_mmcif -.-> cif
  ext_a3m -.-> a3m
  ext_a3m -.-> wol
  ext_rna -.-> rna
  ext_sab -.-> sc40
  ext_sab -.-> sc30
  ext_sig -.-> v1

  ccd --> cif
  cif --> fasta
  cif --> meta
  fasta --> seqid
  seqid --> sc40
  fasta --> sc40
  seqid --> sc30
  fasta --> sc30
  cif --> isi
  seqid --> isi
  isi --> isc
  sc40 --> isc
  fasta --> ppl
  fasta --> tcf
  tcand --> stf

  cif --> att
  seqid --> att
  sc40 --> att
  att --> chain

  a3m --> d16k
  a3m --> d2k
  a3m --> d512

  wol --> hmm
  ppl --> hmm
  hmm --> tcand
  tcf --> tcand
  seqid --> tcand
  fasta --> tcand
  meta --> tcand
  tcand --> tseq
  seqid --> tseq
  fasta --> tseq
  chain --> stm
  tseq --> stm
  stf --> stm
  stm --> ctpl
  tcand --> ctpl

  att --> tr09
  att --> tr26
  tr09 --> mpnn

  cif --> v1
  seqid --> v1
  v1 --> v1a
  seqid --> v1a
  sc30 --> v1a
  v1a --> v2
  sc30 --> v2
  isc --> v2

  classDef build fill:#e9eafb,stroke:#4f5bd5,stroke-width:2px,color:#161a2e;
  classDef rebuild fill:#e0f2f2,stroke:#0e8f8f,stroke-width:2px,color:#0b2b2b;
  classDef materialize fill:#f6ecd8,stroke:#b8791b,stroke-width:2px,color:#2e2408;
  classDef extract fill:#f1e7fb,stroke:#8a4fce,stroke-width:2px,color:#241533;
  classDef parallel fill:#fbe1ea,stroke:#c2416b,stroke-width:2px,color:#33121f;
  classDef external fill:transparent,stroke:#8592a0,stroke-width:1.5px,stroke-dasharray:5 4,color:#6d7a86;
```

## 2. What happens to one node

build-all is itself a datacooker recipe: every node is a step
`name <- submit(*upstream_job_ids)`. The engine topo-sorts (catching cycles / unknown
deps), then threads each node's terminal job id into its dependents.

```mermaid
graph TD
  A["MANIFEST node<br/>(config + upstream job ids)"]:::step
  B["load cfg"]:::step
  C{"output already<br/>on disk?"}:::gate
  S["return None<br/>→ node skipped;<br/>dependents lose one dep"]:::skip
  D{"_op_of(cfg)"}:::gate
  P["submit_pipeline<br/>tiers → merge → index<br/>(internally afterok-chained)"]:::pipe
  W["submit_workflow<br/>single SLURM job"]:::pipe
  J["terminal job_id"]:::done
  DEP["dependents submitted with<br/>--depends-on afterok:job_id"]:::done

  A --> B --> C
  C -->|yes| S
  C -->|no| D
  D -->|build / rebuild| P
  D -->|materialize / extract / parallel| W
  P --> J
  W --> J
  J --> DEP

  classDef step fill:#eef1f4,stroke:#8592a0,stroke-width:1.5px,color:#18222c;
  classDef gate fill:#fff5e6,stroke:#b8791b,stroke-width:2px,color:#2e2408;
  classDef skip fill:#f0f2f4,stroke:#8592a0,stroke-width:1.5px,stroke-dasharray:5 4,color:#4b5a68;
  classDef pipe fill:#e9eafb,stroke:#4f5bd5,stroke-width:2px,color:#161a2e;
  classDef done fill:#e0f2f2,stroke:#0e8f8f,stroke-width:2px,color:#0b2b2b;
```

1. **Skip is native, not bolted on.** The "already built?" gate returns `None`, and
   datacooker's None-propagation drops that node from the run — its dependents simply see
   one fewer upstream job to wait on. No `--condition` flag, no special-casing.
2. **Op decides the executor.** `build`/`rebuild` go to the planning-first pipeline (Ray
   owns memory; tiers → merge → index, each afterok-chained). `materialize`/`extract`/
   `parallel` go to a single workflow job.
3. **SLURM enforces order.** Every node is submitted up front; `--depends-on afterok:<id>`
   makes the scheduler start a node only after its upstreams exit 0. A failing upstream
   cancels its dependents instead of feeding them garbage.

## 3. Op → executor → resources

The op class is inferred from which key the config sets (`_op_of`), and it also fixes the
SLURM footprint for workflow jobs; build/rebuild size themselves from the planning phase.

| Op | Trigger key | Executor | Mem (GB) | Cores | Count |
|----|-------------|----------|---------:|------:|------:|
| `build` | `env_path` | planning-first pipeline | self-sized | self-sized | 6 |
| `rebuild` | `old_env_path` | planning-first pipeline | self-sized | self-sized | 11 |
| `parallel` | `split_recipe_path` | `parallel-run` | 490 | 112 | 2 |
| `extract` | `db_path` + `output_data_path` | `extract-lmdb` | 200 | 32 | 3 |
| `materialize` | `output_data_path` | `workflow run` | 128 | 8 | 9 |

Two `build` nodes (`template/seqid_template_mols`, `template/chain_template`) are `keyed`
— sized by key count rather than byte size.
