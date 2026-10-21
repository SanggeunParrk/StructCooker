#!/bin/bash
#SBATCH -J if_validate
#SBATCH -c 4
#SBATCH --mem=16g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 1:00:00
#SBATCH -o logs/interface/validate_%j.out
# usage: sbatch validate_db.sbatch <db_root> [n] [seed]
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
pixi run python -I scripts/interface/validate_db.py "$@" 2>&1
