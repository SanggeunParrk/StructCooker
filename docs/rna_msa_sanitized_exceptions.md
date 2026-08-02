# RNA MSA Sanitized Exceptions

This records PDB RNA MSA inputs whose RCSB/PDB-derived `polyribonucleotide`
sequence contains non-RNA one-letter symbols. These symbols are present in the
source mmCIF/cif_pdb FASTA, not introduced by the MSA pipeline.

For RNA MSA generation, query sequences are sanitized immediately before
`nhmmer --rna`:

- allowed letters: `A`, `C`, `G`, `U`, `N`, `X`
- all other letters are mapped to `X`

This was needed because `nhmmer --rna` accepts `X` but rejects symbols such as
`F` and `P`.

## Exceptions

| seq_id | original | sanitized | replaced | source PDB chains |
|---|---:|---:|---|---|
| `R00000000000000055022` | `CCAF` | `CCAX` | `F` | `3CPW_EA_.` |
| `R00000000000000024350` | `CCAFX` | `CCAXX` | `F` | `1VQ8_C_.`, `1VQ9_C_.`, `1VQK_C_.` |
| `R00000000000000108229` | `CCAPP` | `CCAXX` | `P` | `5DGV_HF_.`, `5DGV_HF_B`, `5DGV_HF_A`, `5DGV_GF_.`, `5DGV_GF_B`, `5DGV_GF_A` |
| `R00000000000000024349` | `CCAFXX` | `CCAXXX` | `F` | `1VQ6_D_.`, `1VQN_D_.` |
| `R00000000000000016903` | `GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCAF` | `GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCAX` | `F` | `1OB2_B_.` |
| `R00000000000000099359` | `GCCCGGAUAGCUCAGUCGGUAGAGCAGGGGAUUGAAAAUCCCCGUGUCCUUGGUUCGAUUCCGAGUCCGGGCACCAF` | `GCCCGGAUAGCUCAGUCGGUAGAGCAGGGGAUUGAAAAUCCCCGUGUCCUUGGUUCGAUUCCGAGUCCGGGCACCAX` | `F` | `4V5C_Y_.`, `4V5C_DC_.`, `4V5D_V_.`, `4V5D_Y_.`, `4V5D_AC_.`, `4V5D_DC_.`, `4V5J_V_.`, `4V5J_W_.`, `4V5J_CC_.`, `4V5J_DC_.`, `5AFI_Y_A`, `5AFI_Y_B`, `5AFI_Y_.` |
| `R00000000000000016907` | `GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCAFC` | `GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCAXC` | `F` | `1OB5_B_.`, `1OB5_D_.`, `1OB5_F_.` |

## Source Notes

Representative RCSB mmCIF files confirm these are source annotations:

- `3CPW`: entity is `polyribonucleotide` with canonical sequence `CCAF`; residue list is `C C A PHE`.
- `1OB2` and `1OB5`: `TRANSFER-RNA, PHE`, annotated as amino-acylated with `PHE` at the 3' end.
- `1VQ6` and `1VQ8`: synthetic RNA-like polymers include `PHE` and other modified residues whose canonical one-letter sequence includes `F`/`X`.
- `5DGV`: canonical sequence `CCAPP`; `P` is rejected by `nhmmer --rna`, so it is mapped to `X`.

The generated `.a3m` files use the sanitized query sequences.
