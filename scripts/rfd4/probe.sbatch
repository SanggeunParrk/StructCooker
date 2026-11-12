#!/bin/bash
# Conversion probe on a sample of each RFD4-Proteina component. Submit from the repo root:
#   sbatch scripts/rfd4/probe.sbatch [--n 100] [--only tcr_distill ...]
#SBATCH --job-name=rfd4_probe
#SBATCH --partition=cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/rfd4_probe_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
pixi run python scripts/rfd4/probe.py "$@"
