#!/bin/bash
# Exact list of cif_pdb entries with a ligand bonded to two polymer chains. From the repo root:
#   sbatch scripts/unit/bridge_entries_scan.sbatch [OUT_TSV]
#SBATCH --job-name=unit_bridge_scan
#SBATCH --partition=cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH --cpus-per-task=112
#SBATCH --mem=480G
#SBATCH --time=12:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/unit_bridge_scan_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
OUT=${1:-$BIOMOL_ROOT/materials/intermediate/unit/bridge_entries.tsv}
mkdir -p "$(dirname "$OUT")"
# the largest assemblies need a few GB each: 40 workers of a 480 GB node
pixi run python scripts/unit/bridge_entries_scan.py "$OUT" 40
