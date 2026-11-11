#!/bin/bash
# Population numbers and extrapolated sizes for the compact Re-Glyco layouts. Submit from the repo root, after
# reglyco_compact.sbatch has written its JSON and examples.parquet is in materials/intermediate/reglyco_compact/hub:
#   sbatch scripts/analysis/reglyco_extrapolate.sbatch [SAMPLE_JSON]
#SBATCH --job-name=reglyco_extrap
#SBATCH --partition=cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/reglyco_extrap_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
R=$BIOMOL_ROOT/materials/raw/rfd4_proteina/train_datasets/streaming-latest/synthetic/glycans/2026_09_27_reglyco
W=$BIOMOL_ROOT/materials/intermediate/reglyco_compact
SAMPLE=${1:-$W/compact_rs_00005.json}
pixi run python scripts/analysis/reglyco_extrapolate.py --index "$R/blob/index.parquet" --examples "$W/hub/examples.parquet" \
  --sample "$SAMPLE" --out "$W/extrapolation.json"
