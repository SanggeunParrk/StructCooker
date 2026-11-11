#!/bin/bash
# Measure compact layouts for Re-Glyco records on one shipped shard. Submit from the repo root:
#   [LEVELS=6,9,12,19] sbatch scripts/analysis/reglyco_compact.sbatch [SHARD]   (sbatch passes the environment on)
#SBATCH --job-name=reglyco_compact
#SBATCH --partition=cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/reglyco_compact_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
R=$BIOMOL_ROOT/materials/raw/rfd4_proteina/train_datasets/streaming-latest/synthetic/glycans/2026_09_27_reglyco
SHARD=${1:-$R/blob/data/rs_00005.bin}
OUT=$BIOMOL_ROOT/materials/intermediate/reglyco_compact
mkdir -p "$OUT"
pixi run python scripts/analysis/reglyco_compact.py "$SHARD" "$R/blob/index.parquet" --levels "${LEVELS:-6,9,12,19}" --out "$OUT/compact_$(basename "$SHARD" .bin).json"
