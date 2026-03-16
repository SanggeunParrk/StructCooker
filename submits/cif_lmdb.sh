#!/bin/bash
#SBATCH -J cif_DB_급함
#SBATCH --ntasks-per-node=1         
#SBATCH -c 16     
#SBATCH --mem=800g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/cif_lmdb_%a.out
#SBATCH -e ./logs/cif_lmdb_%a.err
#SBATCH --array=0

CONFIG_PATH="configs/cif_lmdb.yaml"
MAP_SIZE=1000000000000 # ~1TB
N_SHARDS=1
export PYTHONPATH="."

python -u scripts/build_lmdb.py build \
  ${CONFIG_PATH} \
  --map-size "${MAP_SIZE}" \
  --shard-idx "${SLURM_ARRAY_TASK_ID}" \
  -n "${N_SHARDS}" \
