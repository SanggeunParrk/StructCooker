#!/bin/bash
#SBATCH -J design_stats
#SBATCH -c 32
#SBATCH --mem=128g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 6:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/design/stats_%j.out
# usage: sbatch build_stats.sbatch --db LMDB --out DIR [--limit N]; then foldseek on DIR/part_*.tar
cd "${SLURM_SUBMIT_DIR:?}"
pixi run python -I scripts/design/build_stats.py "$@" 2>&1 | grep -v WARN
