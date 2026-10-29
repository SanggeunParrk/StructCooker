#!/bin/bash
#SBATCH -J design
#SBATCH -c 4
#SBATCH --mem=32g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 1:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/design/design_%j.out
# usage: sbatch scripts/design/run_py.sbatch scripts/design/check_ss.py --n 100
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
pixi run python -I "$@" 2>&1 | grep -v WARN
