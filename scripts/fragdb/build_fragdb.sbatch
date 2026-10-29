#!/bin/bash
#SBATCH -J frag_build
#SBATCH -c 32
#SBATCH --mem=200g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 6:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/fragdb/fragdb_%j.out
# usage: sbatch build_fragdb.sbatch --out DIR --k 7 --workers 28 --limit 2000
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
pixi run python -I scripts/fragdb/build_fragdb.py "$@" 2>&1 | grep -v WARN
