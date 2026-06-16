#!/bin/bash
#SBATCH -J cif_chain
#SBATCH --ntasks-per-node=1
#SBATCH -c 8
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/%x_%A_%a.out
#SBATCH -e ./logs/%x_%A_%a.err
#SBATCH --array=0

# Usage: N_SHARDS=30 sbatch --array=0-29 submits/maintenance/build_cif_chain.sh
CIF_DB="${CIF_DB:-/data/psk6950/BioMolDB_20260224/cif.lmdb}"
OUT_DB="${OUT_DB:-/data/psk6950/BioMolDB_revisit/cif_chain.lmdb}"
N_SHARDS="${N_SHARDS:-1}"
export PYTHONPATH="."
pixi run python -u -m scripts.maintenance.build_cif_chain \
  --cif-db "${CIF_DB}" --out-db "${OUT_DB}" \
  --shard-idx "${SLURM_ARRAY_TASK_ID}" --n-shards "${N_SHARDS}" --n-jobs 8
