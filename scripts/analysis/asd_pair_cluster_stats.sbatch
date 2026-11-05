#!/bin/bash
#SBATCH -J asd_pstat
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 1
#SBATCH --mem=16g
#SBATCH -t 00:20:00
#SBATCH -o /home/psk6950/StructCooker/logs/staging/asd_pstat_%j.out
# Pair clusters with the antibodies on identical CDRs (seq_cluster_cdr100_ASB.tsv): the antigen as its own cluster
# (seq_id, 100%) and, for comparison, the antigen clusters at 40%. Submit after asd/antibody_seq_cluster_cdr100:
#   sbatch --dependency=afterok:<its job> scripts/analysis/asd_pair_cluster_stats.sbatch
A=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}/asd
S=/home/psk6950/StructCooker/scripts/analysis/asd_pair_cluster_stats.py
echo "=== antigen 100% (seq_id) x antibody CDR-H3/L3 100% ==="
python3 -I $S $A/pairs/pairs_kd1um.tsv seq_id $A/metadata/seq_cluster_cdr100_ASB.tsv $A/pairs/ag100_x_heavy_cdr100_members.tsv
echo
echo "=== antigen 40% x antibody CDR-H3/L3 100% ==="
python3 -I $S $A/pairs/pairs_kd1um.tsv $A/metadata/seq_cluster40_ASD.tsv $A/metadata/seq_cluster_cdr100_ASB.tsv $A/pairs/ag40_x_heavy_cdr100_members.tsv
