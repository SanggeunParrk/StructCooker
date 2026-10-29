#!/bin/bash
# Rebuild everything downstream of the sequence clusters with 40% identity (cluster ids keep their names: cPDB_, cOFD_, cAFM_, cTDM_).
# Run from the repository root. Submits with sbatch dependencies; prints the job ids. Each stage's old output is removed first.
set -euo pipefail
R=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}
SC="pixi run structcooker"
rm -rf $R/metadata/seq_cluster40_{PDB,OFD,AFM,TDM}.tsv $R/metadata/interacting_seq_ids.tsv $R/metadata/interacting_seq_clusters.tsv \
       $R/materials/intermediate/fasta/train.fasta \
       $R/lmdb/pdb/cif/cif_mpnn_train_20210930_res9_chain300.lmdb* $R/lmdb/pdb/valid/cif_attached_valid_1.lmdb* $R/lmdb/pdb/valid/cif_attached_valid_2.lmdb*
# build NAME [deps] -> prints the job id that finishes the stage
build() { local out; out=$($SC build "$1" ${2:+--depends-on "$2"} 2>&1 | grep -v WARN); echo "$out" >&2
          echo "$out" | sed -n 's/^PIPELINE_TERMINAL_JOB=//p; s/^\[slurm\] .*: job //p' | tail -1; }
join() { local IFS=,; echo "$*"; }
J_PDB=$(build metadata/pdb_seq_cluster40)
J_OFD=$(build metadata/ofd_seq_cluster40)
J_AFM=$(build metadata/afm_seq_cluster40)
J_TDM=$(build metadata/teddymer_seq_cluster40)
A_PDB=$(build pdb/cif_attached "$J_PDB")
A_TDM=$(build teddymer/cif_attached "$J_TDM")
A_HOM=$(build afdb_multimer/homodimer_cif_attached "$J_AFM")
A_HET=$(build afdb_multimer/heterodimer_cif_attached "$J_AFM")
A_LONG=$(build distillation/long_cif_attached "$J_OFD")
A_SHORT=$(build distillation/short_cif_attached "$J_OFD")
A_RNA=$(build distillation/rna_cif_attached "$J_OFD")
A_DIS=$(build distillation/disordered_cif_attached "$J_OFD")
D_CON=$(build distillation/disordered_contacts "$A_DIS")
T_21=$(build train/train_20210930 "$A_PDB")
T_26=$(build train/train_20260301 "$A_PDB")
T_FA=$(build metadata/train_fasta "$T_21")
T_MP=$(build train/train_20210930_mpnn "$T_21")
I_IDS=$(build metadata/interacting_seq_ids "$A_PDB")
I_CL=$(build metadata/interacting_seq_clusters "$(join $I_IDS $J_PDB)")
V_AT=$(build valid/valid1_attach "$J_PDB")
V_2=$(build valid/valid2 "$(join $V_AT $I_CL $T_FA)")
echo "clusters: pdb=$J_PDB ofd=$J_OFD afm=$J_AFM tdm=$J_TDM"
echo "attach: pdb=$A_PDB tdm=$A_TDM hom=$A_HOM het=$A_HET long=$A_LONG short=$A_SHORT rna=$A_RNA dis=$A_DIS contacts=$D_CON"
echo "train: 21=$T_21 26=$T_26 fasta=$T_FA mpnn=$T_MP | interacting: ids=$I_IDS clusters=$I_CL | valid: attach=$V_AT valid2=$V_2"
