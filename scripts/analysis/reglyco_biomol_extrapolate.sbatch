#!/bin/bash
# Extrapolate the BioMol-format Re-Glyco layouts to the whole release. Submit from the repo root:
#   sbatch scripts/analysis/reglyco_biomol_extrapolate.sbatch [SHARD_NAME]
#SBATCH --job-name=reglyco_biomol_x
#SBATCH --partition=cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/reglyco_biomol_x_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
R=$BIOMOL_ROOT/materials/raw/rfd4_proteina/train_datasets/streaming-latest/synthetic/glycans/2026_09_27_reglyco
W=$BIOMOL_ROOT/materials/intermediate/reglyco_biomol/${1:-rs_00005}
pixi run python scripts/analysis/reglyco_biomol_extrapolate.py "$W" --index "$R/blob/index.parquet" \
  --population "$BIOMOL_ROOT/materials/intermediate/reglyco_compact/extrapolation.json" --out "$W/extrapolation.json"
