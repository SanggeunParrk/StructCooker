#!/bin/bash
#SBATCH -J pot_frames
#SBATCH -c 32
#SBATCH --mem=200g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 6:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/potential/frames_%j.out
# usage: sbatch count_frames.sbatch <count_frames.py args, e.g. --out DIR --workers 28 --limit 2000>
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
pixi run python -I scripts/potential/count_frames.py "$@" 2>&1 | grep -v WARN
