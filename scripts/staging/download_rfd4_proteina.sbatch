#!/bin/bash
# Parallel download of the selected RFD4-Proteina components: one array task per slice in TASKS.tsv.
#   python3 scripts/staging/download_rfd4_proteina.py manifest      # once, login node (one small API call)
#   N=$(($(wc -l < $DEST/TASKS.tsv) - 1)); sbatch --array=0-$((N-1))%12 scripts/staging/download_rfd4_proteina.sbatch
# Rerun the same command to resume: finished files are skipped, partial ones continue.
#SBATCH --job-name=rfd4_dl
#SBATCH --partition=cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/rfd4_dl_%A_%a.out
set -euo pipefail
# Submit from the repo root: sbatch copies this script, so paths hang off SLURM_SUBMIT_DIR.
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
python3 scripts/staging/download_rfd4_proteina.py fetch --task "$SLURM_ARRAY_TASK_ID" --threads 6
