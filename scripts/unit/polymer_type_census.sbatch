#!/bin/bash
#SBATCH -J unit_poly_census
#SBATCH -c 112
#SBATCH --mem=480g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -t 04:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/analysis/unit_poly_census_%j.out
# usage: sbatch scripts/unit/polymer_type_census.sbatch OUT.tsv      (submit from the repo root)
cd "${SLURM_SUBMIT_DIR:?}"
source scripts/cluster/cssb.env
pixi run python -I scripts/unit/polymer_type_census.py --out "$1" --workers 80 2>&1 | grep -v WARN
