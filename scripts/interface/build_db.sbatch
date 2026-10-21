#!/bin/bash
#SBATCH -J if_build
#SBATCH -c 112
#SBATCH --mem=480g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 12:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/interface/build_%j.out
# usage: sbatch build_db.sbatch <out_dir> <workers> [build_db.py args...]
# multiprocessing needs /dev/shm semaphores: node01/04/07 cannot hold them (see cluster/cssb.env).
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
OUT=${1:?out dir}; W=${2:?workers}; shift 2
pixi run python -I scripts/interface/build_db.py --out "$OUT" --workers "$W" "$@" 2>&1
