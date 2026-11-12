#!/bin/bash
# Check the Re-Glyco delta DB against fresh conversions. Submit from the repo root:
#   sbatch scripts/analysis/reglyco_check.sbatch sample --out DIR [--entries 1]     # dry run on whole entries of one shard
#   sbatch scripts/analysis/reglyco_check.sbatch verify --lmdb $BIOMOL_ROOT/lmdb/glycan/reglyco/reglyco.lmdb --n 5000
#SBATCH --job-name=reglyco_chk
#SBATCH --partition=cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=2-00:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/reglyco_chk_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
pixi run python scripts/analysis/reglyco_check.py "$@"
