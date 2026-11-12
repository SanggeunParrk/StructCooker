#!/bin/bash
# sha256-verify the downloaded RFD4-Proteina files against the manifest (the Hub's LFS digests), one array task per
# download task (TASKS.tsv row): submit from the repo root, after the download finished.
#   N=$(($(wc -l < $BIOMOL_ROOT/materials/raw/rfd4_proteina/TASKS.tsv) - 1)); sbatch --array=0-$((N-1)) scripts/staging/verify_rfd4_proteina.sbatch
#SBATCH --job-name=rfd4_verify
#SBATCH --partition=cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/rfd4_verify_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
T="$BIOMOL_ROOT/materials/raw/rfd4_proteina/TASKS.tsv"
read -r COMP PART < <(awk -F'\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==t {print $2, $3}' "$T")
python3 scripts/staging/download_rfd4_proteina.py verify --sha --component "$COMP" --part "$PART"
