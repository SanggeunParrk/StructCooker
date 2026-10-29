#!/bin/bash
#SBATCH -J tfix_merge
#SBATCH -c 4
#SBATCH --mem=64g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 2:00:00
#SBATCH --exclude=node01,node04,node05,node07
#SBATCH -o logs/template_fix/merge_%j.out
cd "${SLURM_SUBMIT_DIR:?}"
pixi run python -I scripts/template_fix/merge_into_release.py /data/psk6950/template_fix/root /data/shared/cssb_data/BioMol_clean "$@"
