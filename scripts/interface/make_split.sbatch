#!/bin/bash
#SBATCH -J if_split
#SBATCH -c 2
#SBATCH --mem=24g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 1:00:00
#SBATCH -o logs/interface/split_%j.out
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
pixi run python -I scripts/interface/make_split.py "$@" 2>&1
