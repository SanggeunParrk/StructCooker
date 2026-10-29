#!/bin/bash
# The distillation attaches and disordered_contacts, chained to the OFD clustering job.
# usage: rebuild_40_distillation.sh OFD_CLUSTER_JOB
set -euo pipefail
J_OFD=$1
SC="pixi run structcooker"
build() { local out; out=$($SC build "$1" ${2:+--depends-on "$2"} 2>&1 | grep -v WARN); echo "$out" >&2
          echo "$out" | sed -nE 's/^PIPELINE_TERMINAL_JOB=([0-9]+)$/\1/p; s/^\[slurm\] [^:]*: job ([0-9]+).*$/\1/p' | tail -1; }
A_SHORT=$(build distillation/short_cif_attached "$J_OFD")
A_RNA=$(build distillation/rna_cif_attached "$J_OFD")
A_DIS=$(build distillation/disordered_cif_attached "$J_OFD")
D_CON=$(build distillation/disordered_contacts "$A_DIS")
A_LONG=$(build distillation/long_cif_attached "$J_OFD")
echo "attach: short=$A_SHORT rna=$A_RNA dis=$A_DIS contacts=$D_CON long=$A_LONG"
