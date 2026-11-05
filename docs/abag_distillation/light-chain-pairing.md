# Heavy-light pairing of the ASD pairs

Two problems: ASD rows with a **light chain and no heavy chain** (nothing to fold against an antigen as an antibody), and rows with a
**heavy chain and no light chain** that are neither a nanobody nor an scFv (the paratope is VH and VL together; a heavy chain alone folds a binding
site that does not exist, and ipTM does not tell). The second was worked out in the SelfDistillation project (cssb3,
`/public_data/kkh517/SelfDistillation/docs/Light-chain-pairing.md`, `Paired-set-build.md`); this page is what we took from it, what we measured, and
what we do.

## What SelfDistillation did

ASD has no antibody identifier: a patent's 50 antibodies stored as half rows lose their pairing for good (heavy only 429,902 rows, light only 79,105,
the other chain's numbering is empty too: the chain was never there upstream). Its fix, rules in order, each row keeping `chain_status` and its source:

1. `paired_original` (ASD has the light chain), 2. `scfv` (a `(G{3,5}S){2,}` linker in the heavy chain, checked against the flag: 37 rows had the flag and no
linker), 3. `vhh` (nanobody flag), 4. `filled_covid` (AntiBinder CoV-AbDab: exact (heavy, antigen) match, one light chain), 5. `filled_hiv` (the heavy chain
is a prefix of AntiBinder's `antibody`, the rest >= 60 aa is the light chain), 6. `filled_biomap` (the pair column), 7. `filled_db` (POISK + PLAbDab: exact heavy chain, **exactly
one** light-chain candidate), 8. `ambiguous` (two or more: not filled), 9. `unresolved`. The principle: a wrong light chain is worse than none.
Their result: the 132,747 pairs without a light chain became 61,428 on hold (54% fewer).

## What we measured (the ASD pairs, 535,295 of them at the time; SAbDab as ground truth)

| | pairs |
|---|---:|
| light chain in ASD (`paired_original`) | 272,423 |
| scFv | 78,679 |
| VHH | 91,822 |
| heavy chain only | 92,371 |

SAbDab (8,010 PDB entries, 5,184 distinct heavy-light pairs taken from the crystal structures) is a ground truth for pairing. A pairing counts as
consistent when its light chain is, or contains, a light chain SAbDab pairs with that heavy chain (containment either way: the PDB chain carries a
constant domain or a tag). `scripts/abag/asd_pairing_gold.py`:

| pairing | heavy chain in SAbDab | consistent |
|---|---:|---:|
| ASD `structures-antibodies` (control: these are PDB entries) | 1,518 | 99.8% |
| ASD `literature` | 1,617 | 97.3% |
| AntiBinder covid-19 | 382 | 94.8% |
| ASD `patents` | 1,828 | 86.8% |
| POISK, heavy chain with exactly one light chain | 629 | 86.2% |
| heavy-only pairs filled from POISK (unique candidate) | 7,248 | 85.6% |

- **The `filled_db` rule is wrong about 14% of the time, not 3%.** SelfDistillation's 2.9% is the disagreement of POISK and PLAbDab on heavy chains
  both name; both are made from patents, so it is not an independent check. Against crystal structures a unique POISK candidate disagrees in 14%
  (the same as ASD's own patent pairing, 13%). Part of that is a legitimate other pairing (a patent lists several light-chain variants of one heavy chain),
  but for choosing the light chain of *this* antibody it is an error all the same.
- **Requiring the antigen too does not help.** POISK carries the target sequence for only about half of its rows; among the 2,742 pairs where a POISK row with
  the same heavy chain and antigen exists, ASD's light chain is among its candidates in 23%, and the unique candidate agrees in 87% (519 of 594).
- **AntiBinder is the sturdy source.** covid-19 and hiv are the same rows ASD took, so the join is exact (SelfDistillation: 26,864 of 26,864 (heavy, antigen)
  pairs; hiv heavy chain a prefix of the original in 190 of 192; antigen 940 of 940), and the covid-19 light chains agree with SAbDab in 94.8% (362 of 382).
  biomap is weaker: 68% (408 of 598), a PDB-derived set with full-chain Fabs; it is 896 pairs.
- The control and the variant libraries: `flab_koenig2017` (2,014 pairs) reads 0% because its light chains are mutants and the test asks for the wild-type
  chain; the test is for natural pairings, not for libraries.

## What we do

Stage `abag/pairs_chains` (recipe `workflows/metadata/abag_pairs_chains.py`, code `instructions/readers/abag_chains.py`) adds `chain_status`, `light_source` and, for a recovered chain, `light_sequence` to the
pair table (`abag/pairs/pairs_chains.tsv`), for the pairs of all three sources. The rules are SelfDistillation's 1-6, 8, 9; **POISK does not fill** (rule 7 is left out; PLAbDab is not here and has no
antigen). The flags are the `nanobody` / `scfv` columns of the pair table (ASD's flags, POISK's `VHH` type); AbRank has no heavy-only pairs.

| `chain_status` | pairs (677,817) |
|---|---:|
| `paired_original` | 412,781 (after merging, in `pairs_basic`) |
| `vhh` | 92,331 |
| `scfv` | 78,680 |
| `filled_covid` | 13,040 |
| `filled_hiv` | 12,040 |
| `filled_biomap` | 921 |
| `unresolved` | 66,431 (ASD `patents` 62,693, `genbank` 1,717, `met` 1,516, POISK 358, ...) |
| `ambiguous` | 293 |

**Default selection: keep everything except `unresolved` and `ambiguous` (66,724 pairs out).** Light-only rows are already gone (S2 needs a heavy chain: 66,668 ASD rows, 66,231 of them `patents`).
About 22,000 of the unresolved ASD `patents` pairs have a unique POISK candidate; they are out by default (about one in seven would carry the wrong light chain) and can be brought back with a
later POISK stage.

A recovered light chain is a new sequence: `light_seq_id` is empty until the final seq_id stage (`abag/final_seq_id_map`) takes `light_sequence` up. ASD's own patent pairings stay (13%
inconsistent by the control above): drop `asd:patents` pairs without a second source if that is too much.

## Not done

PLAbDab (176,894 pairs, no antigen, licence unrecorded); the patent full texts (about 38,000 patents, an LLM extraction); folding the filled pairs with
two chains and with three and comparing ipTM, which is what would say whether the recovery matters for the structures.
