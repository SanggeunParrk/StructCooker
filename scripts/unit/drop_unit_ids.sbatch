#!/bin/bash
# Remove stale unit ids after a unit_pdb patch. From the repo root:
#   sbatch --dependency=afterok:<patch apply job> scripts/unit/drop_unit_ids.sbatch IDS.txt
#SBATCH --job-name=unit_drop_ids
#SBATCH --partition=cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH --output=logs/unit_drop_ids_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/cluster/cssb.env
pixi run python scripts/unit/drop_unit_ids.py "$BIOMOL_ROOT/lmdb/unit/pdb/unit_pdb.lmdb" "$1"
