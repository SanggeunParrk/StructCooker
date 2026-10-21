#!/bin/bash
#SBATCH -J if_examples
#SBATCH -c 2
#SBATCH --mem=24g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 0:30:00
#SBATCH -o logs/interface/examples_%j.out
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
pixi run python -I scripts/interface/pick_examples.py "$@" 2>&1
