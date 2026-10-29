#!/bin/bash
# valid2 and the distillation attaches, after rebuild_40_rest2.sh submitted the rest.
# usage: rebuild_40_rest3.sh VALID1_ATTACH_JOB INTERACTING_CLUSTERS_JOB TRAIN_FASTA_JOB OFD_CLUSTER_JOB
set -euo pipefail
V_AT=$1; I_CL=$2; T_FA=$3; J_OFD=$4
SC="pixi run structcooker"
build() { local out; out=$($SC build "$1" ${2:+--depends-on "$2"} 2>&1 | grep -v WARN); echo "$out" >&2
          echo "$out" | sed -nE 's/^PIPELINE_TERMINAL_JOB=([0-9]+)$/\1/p; s/^\[slurm\] [^:]*: job ([0-9]+).*$/\1/p' | tail -1; }
V_2=$(build valid/valid2 "$V_AT,$I_CL,$T_FA")
A_LONG=$(build distillation/long_cif_attached "$J_OFD")
A_SHORT=$(build distillation/short_cif_attached "$J_OFD")
A_RNA=$(build distillation/rna_cif_attached "$J_OFD")
A_DIS=$(build distillation/disordered_cif_attached "$J_OFD")
D_CON=$(build distillation/disordered_contacts "$A_DIS")
echo "valid2=$V_2 | attach: long=$A_LONG short=$A_SHORT rna=$A_RNA dis=$A_DIS contacts=$D_CON"
