#!/bin/bash
#SBATCH -J design_db
#SBATCH -c 56
#SBATCH --mem=300g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 24:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/design/db_%j.out
# usage: sbatch build_db.sbatch --db SRC --items items.tsv --out DIR --source NAME [--part i/n] [--limit N]
cd "${SLURM_SUBMIT_DIR:?}"
pixi run python -I scripts/design/build_db.py "$@" 2>&1 | grep -v WARN
