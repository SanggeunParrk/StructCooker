#!/bin/bash
#SBATCH -J pdb_special
#SBATCH -c 112
#SBATCH --mem=480g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -t 24:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/analysis/pdb_special_%j.out
# usage: sbatch pdb_special_chains.sbatch OUT [--ids a,b,c] [--ccd-only] [--every K] [--max-pep-len N]
cd "${SLURM_SUBMIT_DIR:?}"
OUT=$1; shift
pixi run python -I scripts/analysis/pdb_special_chains.py --out "$OUT" --workers 80 "$@" 2>&1 | grep -v WARN
