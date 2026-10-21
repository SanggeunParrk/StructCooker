#!/bin/bash
#SBATCH -J if_finalize
#SBATCH -c 2
#SBATCH --mem=48g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 3:00:00
#SBATCH -o logs/interface/finalize_%j.out
# usage: sbatch finalize.sbatch <db_root>
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
pixi run python -I scripts/interface/finalize.py "$@" 2>&1
