#!/bin/bash
# Key-set check: every attached DB holds exactly the keys of its source (or a subset, where entries are filtered). sbatch only.
L=/data/shared/cssb_data/BioMol_clean/lmdb
keys() { grep -v -e '^#' -e '^key	' "$1.index.tsv" | cut -f1 | LC_ALL=C sort; }
cmp() {  # name source attached
  keys $L/$2 > /tmp/ks_a.$$; keys $L/$3 > /tmp/ks_b.$$
  a=$(wc -l < /tmp/ks_a.$$); b=$(wc -l < /tmp/ks_b.$$)
  only_src=$(LC_ALL=C comm -23 /tmp/ks_a.$$ /tmp/ks_b.$$ | wc -l); only_att=$(LC_ALL=C comm -13 /tmp/ks_a.$$ /tmp/ks_b.$$ | wc -l)
  echo "$1: source $a attached $b | only in source $only_src | only in attached $only_att"
}
cmp afm_hetero afdb_multimer/cif/cif_afm_heterodimer.lmdb afdb_multimer/cif_attached/cif_afm_heterodimer_attached.lmdb
cmp afm_homo afdb_multimer/cif/cif_afm_homodimer.lmdb afdb_multimer/cif_attached/cif_afm_homodimer_attached.lmdb
cmp teddymer teddymer/cif/cif_teddymer.lmdb teddymer/cif_attached/cif_teddymer_attached.lmdb
cmp short distillation_short/cif/cif_short.lmdb distillation_short/cif_attached/cif_short_attached.lmdb
cmp rna distillation_rna/cif/cif_rna.lmdb distillation_rna/cif_attached/cif_rna_attached.lmdb
cmp disordered distillation_disordered/cif/cif_disordered.lmdb distillation_disordered/cif_attached/cif_disordered_attached.lmdb
cmp pdb pdb/cif/cif_pdb.lmdb pdb/cif_attached/cif_pdb_attached.lmdb
cmp valid_1 pdb/valid/cif_valid_1.lmdb pdb/valid/cif_attached_valid_1.lmdb
cmp pdb_train_in_pdb pdb/cif_attached/cif_pdb_attached.lmdb pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb
cmp valid2_in_valid1 pdb/valid/cif_attached_valid_1.lmdb pdb/valid/cif_attached_valid_2.lmdb
cmp long distillation_long/cif/cif_long.lmdb distillation_long/cif_attached/cif_long_attached.lmdb
