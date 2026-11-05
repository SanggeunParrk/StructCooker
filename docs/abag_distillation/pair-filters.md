# Basic filters of the pairs

The pair table is built in stages (each a `structcooker build abag/<name>`), every filter a column or a documented rule, so nothing disappears without a trace.
Order and counts (ASD + AbRank + POISK; a pair is one (heavy, light, antigen) sequence triple):

| # | Filter | Stage | Pairs left |
|---|---|---|---:|
| | S2 + antigen 10-3000 aa + S3 + Kd <= 1 uM, three sources (no S1: structures are handled in step 4) | `abag/pairs_binding` | 677,817 |
| 1 | duplicates merged | `abag/pairs_basic` | 676,517 |
| 2 | usable chain status (light-only rows are gone already) | `abag/pairs_chains`, `abag/pairs_basic` | 609,793 |
| 3 | polyprotein | `abag/pairs_basic` | 609,573 |
| 4 | an experimental structure with the same heavy chain, light chain and antigen | `abag/pairs_basic` | **607,342** |
| 5 | the antigen has an MSA of depth over 50 or a template (2021-09-30) | `abag/ag_keep`, `abag/final_*` | see [RECIPE.md](RECIPE.md) |

Output `$BIOMOL_ROOT/abag/pairs/pairs_basic.tsv`, one row per pair with all columns (`antigen_seq_id, heavy_seq_id, light_seq_id, light_sequence, chain_status, light_source, datasets,
polyprotein, pdb_structure, keep`); `keep == 1` is the set after steps 1-4. The filters are columns, so a different choice is a different `keep`, not a rebuild.
`datasets` is `source:dataset` (`asd:patents`, `abrank:AlphaSeq`, `poisk:Conventional`). By source after step 4: ASD 484,207, AbRank 122,984, POISK 10,593 (473,899 only ASD, 113,325 only AbRank,
10,175 only POISK, about 10,000 in two sources).

## Before the stage: what is already out

- **S2**: a row needs a heavy chain and an antigen, all chains A-Z. **Light-only rows are out here** (ASD: 66,668 rows, 66,231 of them `patents`).
- **Antigen length**: 10-3000 aa (this is where polyproteins over 3,000 aa and the giant single chains leave; peptides under 10 aa have nothing to align). AbRank antigens of several chains
  (`:`) are left out.
- **Kd <= 1 uM**: [affinity-filter.md](affinity-filter.md), per source.
- **S1 is off**: ASD's `target_pdb` is the structure of the *parent* (the Buzz library's 524,346 rows all carry the trastuzumab-HER2 PDB id), so it cannot say whether this antibody has a structure. Step 4
  asks that directly.

## 1. Duplicates

A pair is one sequence triple, so a row repeated in a source, or the same pair in two datasets or sources, is one pair (its `datasets` lists both): that is `abag/pairs_binding`.
`abag/pairs_basic` merges once more, because filling a light chain (`abag/pairs_chains`) can produce a triple that a source also holds: **1,300 pairs merged**, datasets joined,
and a recovered light chain that is a sequence we already hold takes that seq_id.

Not done: near-duplicates (a point mutant is another pair; the alphaseq and Buzz libraries are variant libraries): the clusters and a per-cluster cap are the way to deal with them, in a later step.
The antigens differ between sources (ASD holds the UniProt chain, AbRank the construct, POISK the patent's fragment), so the same antibody-target pair from two sources is mostly two pairs.

## 2. Chain status

[light-chain-pairing.md](light-chain-pairing.md). `keep` takes `paired_original`, `scfv`, `vhh` and the three `filled_*`; `unresolved` and `ambiguous` (66,724 pairs: ASD 66,388, POISK 358) are out.

## 3. Polyprotein

An antigen is a polyprotein when a target name of it (ASD `target_name`, AbRank `Ag_name`, POISK `target.protein_name`) matches `polyprotein`, `ORF1a/ORF1ab`, `pp1a/pp1ab` or `gag-pol`
(`POLYPROTEIN` in `abag_basic.py`): **220 pairs**. Why: a polyprotein is cut into products that fold on their own, so the uncut chain is not one folding unit and the antibody binds a product;
SelfDistillation planned to bring them back as their cleavage products (not done here). The length rule is not a second filter: an antigen over 3,000 aa is not in the table. The 800 aa limit of
SelfDistillation's later funnel is not applied.

## 4. Experimental structure

Remove a pair when a PDB complex has **the same heavy chain, the same light chain and the same antigen, all identical (100% sequence identity: the same seq_id for every chain)**.
A structure of the **parent does not remove the variant**: the Buzz and alphaseq libraries are point mutants of crystallised antibodies and have none of their own, and ASD's
`target_pdb` (set on 247,046 of the ASD pairs, nearly all of it a parent's structure) is not used.

How (`PdbComplexIndex.find(..., match="exact")`): SAbDab (`materials/raw/sabdab/sabdab_summary_all.tsv`, entries with a heavy chain and a protein or peptide antigen) gives the chains of each
entry, their sequences come from the PDB fasta (`cif_pdb.fasta`). A pair has a structure when one entry has its heavy chain, its light chain (none, for a single-chain pair: the entry then
has no light chain either) and one of its antigen chains, each sequence equal. The PDB id goes in the `pdb_structure` column.

Result: **2,231 pairs** flagged (asd:structures-antibodies 1,526, asd:structures-nanobodies 583, asd:biomap 353, asd:skempiv2 12, asd:literature 7, abrank:SabDab 4, abrank:SKEMPIv2 3,
asd:ab-bind 3). **Control**: the ASD datasets built from PDB entries (`structures-antibodies`, `structures-nanobodies`) have 2,146 pairs, of which 2,107 (98.2%) are flagged. Buzz (172,149 pairs) and
alphaseq (128,973) are not flagged at all, as they should be.

What exact identity leaves in. The looser reading (`structure_match: contained`: an antibody chain contained in a PDB chain or the other way, an antigen that is a fragment of the
PDB's or has 80% of its 10-mers in it) flags about twice as many pairs (4,616 against 2,204 on the ASD pairs alone). The extra ones are the same antibody in another construct: a variable
domain where the PDB holds a Fab chain, the full UniProt chain where the PDB holds a fragment. Under the exact rule they stay in the set, although a model has seen the same interface in
training. It is one parameter (`structure_match` in `db_configs/abag/pairs_basic.yaml`).

Limits: only SAbDab entries count (a structure SAbDab does not list is missed; the summary is the 2026-03-16 one); the test is on sequences, not on the epitope, and an exact match needs the
residues of the construct to agree (an expression tag or a modified residue in the PDB chain makes it a miss).

## 5. The antigen filter

After the MSA and the templates of the antigens are searched ([sources.md](sources.md)): the antigen is dropped when its MSA depth is 50 or less and it has no template in the 2021-09-30 version
(`abag/worklists/ag_keep.txt`). The final FASTAs, seq_ids, clusters and the final pair table only use the kept pairs of kept antigens (`run.sh final`, [RECIPE.md](RECIPE.md) step 4).
