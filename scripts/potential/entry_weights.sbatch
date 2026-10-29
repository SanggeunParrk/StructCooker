#!/bin/bash
#SBATCH -J pot_weights
#SBATCH -c 2
#SBATCH --mem=48g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 1:00:00
#SBATCH -o logs/potential/weights_%j.out
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
pixi run python -I scripts/potential/entry_weights.py "$@" 2>&1 | grep -v WARN
