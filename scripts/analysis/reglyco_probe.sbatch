#!/bin/bash
# Compression probe on one downloaded Re-Glyco shard. Submit from the repo root:
#   sbatch scripts/analysis/reglyco_probe.sbatch [SHARD]
#SBATCH --job-name=reglyco_cmp
#SBATCH --partition=cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/reglyco_probe_cmp_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
SHARD=${1:-$BIOMOL_ROOT/materials/raw/rfd4_proteina/train_datasets/streaming-latest/synthetic/glycans/2026_09_27_reglyco/blob/data/rs_00005.bin}
pixi run python scripts/analysis/reglyco_probe.py "$SHARD"
