#!/bin/bash
# Scan every unit of lmdb/unit/pdb for polymer nodes that hold more than one chain. From the repo root:
#   sbatch scripts/unit/merged_node_scan.sbatch [OUT_TSV]
#SBATCH --job-name=unit_merged_scan
#SBATCH --partition=cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/unit_merged_scan_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
OUT=${1:-$BIOMOL_ROOT/materials/intermediate/unit/merged_nodes.tsv}
mkdir -p "$(dirname "$OUT")"
pixi run python scripts/unit/merged_node_scan.py "$BIOMOL_ROOT/lmdb/unit/pdb" "$OUT" 44
