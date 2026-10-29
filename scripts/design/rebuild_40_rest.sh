#!/bin/bash
# Second half of rebuild_40.sh, after pdb/cif_attached was resubmitted alone (it shares logs/pipeline/cif_attached with teddymer/cif_attached,
# so the two must never be submitted together). usage: rebuild_40_rest.sh PDB_ATTACH_JOB OFD_CLUSTER_JOB
set -euo pipefail
A_PDB="${1:-}"; J_OFD=$2
SC="pixi run structcooker"
build() { local out; out=$($SC build "$1" ${2:+--depends-on "$2"} 2>&1 | grep -v WARN); echo "$out" >&2
          echo "$out" | sed -n 's/^PIPELINE_TERMINAL_JOB=//p; s/^\[slurm\] .*: job //p' | tail -1; }
join() { local IFS=,; echo "$*"; }
T_21=$(build train/train_20210930 "$A_PDB")
T_26=$(build train/train_20260301 "$A_PDB")
T_FA=$(build metadata/train_fasta "$T_21")
T_MP=$(build train/train_20210930_mpnn "$T_21")
I_IDS=$(build metadata/interacting_seq_ids "$A_PDB")
I_CL=$(build metadata/interacting_seq_clusters "$I_IDS")
V_AT=$(build valid/valid1_attach)
V_2=$(build valid/valid2 "$(join $V_AT $I_CL $T_FA)")
A_LONG=$(build distillation/long_cif_attached "$J_OFD")
A_SHORT=$(build distillation/short_cif_attached "$J_OFD")
A_RNA=$(build distillation/rna_cif_attached "$J_OFD")
A_DIS=$(build distillation/disordered_cif_attached "$J_OFD")
D_CON=$(build distillation/disordered_contacts "$A_DIS")
echo "attach: long=$A_LONG short=$A_SHORT rna=$A_RNA dis=$A_DIS contacts=$D_CON"
echo "train: 21=$T_21 26=$T_26 fasta=$T_FA mpnn=$T_MP | interacting: ids=$I_IDS clusters=$I_CL | valid: attach=$V_AT valid2=$V_2"
