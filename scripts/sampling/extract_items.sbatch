#!/bin/bash
#SBATCH -J extract_items
#SBATCH -c 112
#SBATCH --mem=400g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -t 3-00:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/design/extract_%j.out
# usage: sbatch extract_items.sbatch DB OUT [--part i/n]   (cpu-long: preemptible, hence the clean start)
cd "${SLURM_SUBMIT_DIR:?}"
DB=$1; OUT=$2; shift 2
rm -rf "$OUT"   # a requeued (preempted) job starts clean: the workers append
pixi run python -I scripts/sampling/extract_items.py --db "$DB" --out "$OUT" --workers 112 "$@" 2>&1 | grep -v WARN
