# Binding filter: Kd <= 1 uM

One rule for every source: **the antibody binds when its Kd is at most 1 uM.** The sources do not report Kd. ASD, AbRank and
POISK carry binary labels, assay scores, neutralisation concentrations, free energies, ratios and free text, so each type is
read on its own and turned into one verdict per row: binds / does not bind / **no verdict**. A pair (heavy, light, antigen)
is kept when some row says it binds and no row says it does not. A pair whose rows all give no verdict is dropped.

Code: `structcooker/instructions/readers/asd_binding.py` (`binding_verdict`, `pair_binds`), used by the stage
`abag/pairs_binding` ([sources.md](sources.md)) for all three sources: ASD in `asd_binding.py`, AbRank and POISK in `abag_binding.py`. Pairs by source after the rule: ASD 552,930, AbRank 122,990, POISK 10,822 (a pair two sources report counts in both).

Rows are counted in ASD v2.0 (1,177,819 rows) before any filter.

## ASD v2.0 (implemented)

| `affinity_type` | datasets (rows) | binds when | notes |
|---|---|---|---|
| `-log Kd (M)` | abbd (152,401), flab_koenig2017 (4,275) | pKd >= 6 | abbd 5.6-12.4, koenig 7.2-12.2 |
| `Kd (nM)` | htm (361), absci (354), aatp (93), aae (35) | <= 1000 nM | |
| `Kd (M)` / `KD (M)` | scmc (754), seha (15) | <= 1e-6 M | scmc: 462 `nan` and 102 `I.C.` have no verdict, 190 numeric |
| `delta_g` | biomap (2,725) | <= -8.185 kcal/mol | **assumed kcal/mol and 25 C (298.15 K)**: dG = RT ln(Kd), RT = 0.5925 kcal/mol; the release states neither |
| `EC50 (nM)` | htm (21) | <= 1000 nM | a half-effect concentration is about a Kd or above it, so <= 1 uM says Kd <= ~1 uM; above 1 uM gives no verdict (all 21 are 0.63-303 nM) |
| `alphaseq` | alphaseq (198,703) | score <= 3.0 | log10(Kd in nM) on the AlphaSeq scale, lower binds tighter, 3.0 = 1 uM; 8,229 `nan` have no verdict. The scale is the public `aalphabio/open-alphaseq` one; that ASD copied it unchanged is inferred from the range (-1.0 .. 6.7), not checked |
| `bool` | patents, hiv, covid-19, literature, genbank, met, absci, structures-*, ... | label 1 / `True` | the label, not a Kd: **1 uM is not checked**. `0` / `False` is a stated non-binder and drops the pair. patents, literature, genbank, structures-* are all 1: positives only |
| `fuzzy` | buzz (524,346) | `h` | FACS bins of the Buzz trastuzumab-HER2 library (Chinery et al.): the high bin binds, medium and low do not ("we assign 'high' ... positive binders", "group 'medium' and 'low' ... negative binders"). The paper never writes h / m / l; h = high is inferred from 172,149 / 524,346 = 32.8%, its stated share of high-affinity binders. No Kd behind it |

### Left out as ambiguous: no verdict, so these rows never make a pair

| `affinity_type` | dataset (rows) | why |
|---|---|---|
| `-log Kd (M)` | **flab_warszawski2019 (2,048)** | values -1.78 .. 4.43. Read as pKd in molar they are millimolar non-binders; read as -log10(Kd in nM) they would be pM binders. The label fits neither reading, so the dataset is excluded by name (`UNIT_UNCERTAIN_DATASETS`). Before this rule all 2,048 rows were judged non-binders |
| `-log Kd (nM)` | flab_shanehsazzadeh2023 (448), flab_hie2023 (146), flab_rosace2023 (19) | the label says nM but the values of the first two (5.8-9.0, 4.5-9.5) look like pKd in molar; psh's cross-check on cssb-master (AbRank 31 nM against ASD 7.51 for the same measurement) supports that for flab_shanehsazzadeh2023 only. rosace2023 (-0.38 .. 1.59) fits neither. Not judged on a guess |
| `kd` | osh (30) | no unit |
| `ka (1/Ms)` | htm (3) | an on-rate, not an affinity |
| `ic_50` | dlgo (360) | a neutralisation / inhibition concentration with **no unit** in the release (0.001-5.7, 11 rows `>8.8`), so the ug/mL rule below cannot be applied |
| `ddg`, `log_enrichment`, `elisa_mut_to_wt_ratio` | skempiv2 (400), ab-bind (270), abbd (3,452), abdesign (658) | relative to a parent: no absolute Kd |

The tables above are the decisions taken so far: **ambiguous types are dropped, not guessed.** To bring one back, give its
unit evidence and add the rule to `binding_verdict`.

## AbRank (implemented)

`AbRank_dataset_updated_v2.csv`, 342,362 rows; the unit is in the column header. The rules `abag_binding.abrank_verdict` applies.

| column | rows | what it is | rule |
|---|---:|---|---|
| `Affinity_Kd [nM]` | 94,619 | Kd | <= 1000 nM binds, above does not. 94,015 (99.4%) are already <= 1 uM, so the rule barely filters |
| `IC50 [ug/mL]` with `Aff_op` `=` | 40,354 (CATNAP 39,038, AbCoV 1,316) | neutralisation concentration | **<= 10 ug/mL binds; above gives no verdict** (89% of CATNAP's 39,038 pass) |
| `Aff_op` `>` / `<` (CATNAP 33,984) | | IC50 beyond the tested range | no verdict: not neutralising is not "does not bind" |
| `escape`, `log(Kd_ratio)` | RBD-escape 192,559 | loss of binding of an RBD mutant, 0-3.6 and 0-2 | **no verdict**, except the 20,438 rows that also have `Affinity_Kd [nM]`, which follow the Kd rule |

Why these:

- **IC50 <= 10 ug/mL.** Neutralising needs binding, and an IC50 is about the Kd or above it (it cannot sit far below it), so a small IC50 is
  sufficient evidence of a small Kd. 1 uM is 150 ug/mL for an IgG (150 kDa), 50 for a Fab, 27 for an scFv, 15 for a VHH; **10 ug/mL is below
  all of them**, so one number is safe for every format. A larger IC50 is weak neutralisation, not proof of non-binding: no verdict, not "no".
  <= 150 would keep 99.7% but is right only for an IgG; <= 1 would keep 65% and discards good binders.
- **Escape gives no verdict.** The escape fraction (Bloom-lab deep mutational scanning) says how much of an RBD mutant's
  yeast population lost antibody binding at one antibody concentration; it says that the mutation matters, not what the Kd is. The only rows
  that could calibrate a threshold are the 20,438 that also have a Kd, and 20,433 of them are <= 1 uM: no negatives, so no cut could
  be tested. `escape <= 0.05` (115,951 rows, 60%) would have been a guess.

AbRank's `AlphaSeq` source is a Kd in nM (71,834 rows, median 3.4 nM, max 10,275), not the ASD `alphaseq` score.

## POISK (implemented)

`target.binding_affinity.affinity_value` is free text; `audit.verdict` is POISK's own check against the patent. 45,435 of 89,074 rows have a heavy chain and a target sequence.

- 42,614 rows have a value; 35,899 read as a number with a molar unit. Of those, 1,619 are labelled Kd (909 <= 1 uM), 1,986 IC50, 1,331 EC50, 69 Ki, and about 30,900 carry **no label** (21,037 <= 1 uM).
- Texts that are not a Kd: `Titer 1:64`, `Not Detectable`, `4-fold lower than chimeric IgG1`, `8.77x10^9 M-1` (an association constant), `10^-7 to 10^-8 M` (a range).
- **Rules** (`poisk_verdict`):
  - `audit.verdict` `CONTRADICTED` or `NOT_FOUND`: the row is dropped (14,665 + 91 rows).
  - **No affinity text: the patent reports the pair as binding, a binder** -- the same status as the ASD `patents` rows (`bool = 1`), which come from the same kind of source.
  - A negative phrase (`Not Detectable`, `no binding`): does not bind.
  - A value **labelled Kd** in a molar unit: <= 1 uM binds, above does not (`<` and `>` count only when they decide it).
  - A measured `IC50` / `EC50`: molar <= 1 uM or ug/mL <= 10 binds, above gives no verdict.
  - An unlabelled concentration, `Ki`, `M-1`, a range, a ratio to another antibody, a titer: **no verdict**: a bare "5 nM" in a patent can be a Kd, an IC50 or a titration point.
- The earlier proposal gave a row without text no verdict; that dropped the largest group of POISK pairs while the same kind of ASD row was kept, so it was changed.

## Assumptions to keep in mind

- **Positive-only datasets** (patents 190,484, literature, genbank, structures-*) cannot be checked against a Kd; their label is
  taken as binding.
- **Mixed signals.** A pair reported twice with a binder row and a non-binder row is dropped ("some row yes, no row no").
- **A binder is not a measured Kd <= 1 uM** for `bool` and `fuzzy` rows: it is the dataset's own label.
