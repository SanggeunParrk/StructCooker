#!/bin/bash
# Per-residue ig_region / cdr table of an RFD4-Proteina set. Submit from the repo root:
#   sbatch scripts/rfd4/residue_annotations.sbatch antibody_distill|tcr_distill
#SBATCH --job-name=rfd4_resid
#SBATCH --partition=cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/rfd4_resid_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
pixi run python scripts/rfd4/residue_annotations.py "$1" --workers "$SLURM_CPUS_PER_TASK"
