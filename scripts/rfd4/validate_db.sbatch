#!/bin/bash
# Validate a built RFD4-Proteina set DB. Submit from the repo root:
#   sbatch scripts/rfd4/validate_db.sbatch <set> [--sample 20000]
#SBATCH --job-name=rfd4_validate
#SBATCH --partition=cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/rfd4_validate_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
pixi run python scripts/rfd4/validate_db.py "$@"
