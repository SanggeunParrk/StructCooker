#!/bin/bash
#SBATCH -J msa_remap
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=80g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -w node02
#SBATCH -o ./logs/msa_remap.out
#SBATCH -e ./logs/msa_remap.err

set -euo pipefail

SEQ_ID_MAP_PATH="$1"
OLD_MSA_DIR="$2"
NEW_MSA_DIR="$3"

python -u scripts/maintenance/remap_msa.py   "$SEQ_ID_MAP_PATH"   "$OLD_MSA_DIR"   "$NEW_MSA_DIR"
