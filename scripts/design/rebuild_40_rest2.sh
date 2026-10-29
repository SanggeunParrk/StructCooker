#!/bin/bash
# Stages that can be planned now (their inputs exist). train_20210930_mpnn reads train_20210930, which must be built first, so it is submitted afterwards.
# usage: rebuild_40_rest2.sh TRAIN_FASTA_JOB OFD_CLUSTER_JOB
set -euo pipefail
T_FA=$1; J_OFD=$2
SC="pixi run structcooker"
build() { local out; out=$($SC build "$1" ${2:+--depends-on "$2"} 2>&1 | grep -v WARN); echo "$out" >&2
          echo "$out" | sed -nE 's/^PIPELINE_TERMINAL_JOB=([0-9]+)$/\1/p; s/^\[slurm\] [^:]*: job ([0-9]+).*$/\1/p' | tail -1; }
join() { local IFS=,; echo "$*"; }
I_IDS=$(build metadata/interacting_seq_ids)
I_CL=$(build metadata/interacting_seq_clusters "$I_IDS")
V_AT=$(build valid/valid1_attach)
V_2=$(build valid/valid2 "$(join $V_AT $I_CL $T_FA)")
A_LONG=$(build distillation/long_cif_attached "$J_OFD")
A_SHORT=$(build distillation/short_cif_attached "$J_OFD")
A_RNA=$(build distillation/rna_cif_attached "$J_OFD")
A_DIS=$(build distillation/disordered_cif_attached "$J_OFD")
D_CON=$(build distillation/disordered_contacts "$A_DIS")
echo "interacting: ids=$I_IDS clusters=$I_CL | valid: attach=$V_AT valid2=$V_2"
echo "attach: long=$A_LONG short=$A_SHORT rna=$A_RNA dis=$A_DIS contacts=$D_CON"
