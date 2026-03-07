#!/bin/bash
#SBATCH -J ccd_DB
#SBATCH --ntasks-per-node=1         
#SBATCH -c 80       
#SBATCH --mem=400g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/ccd_lmdb_%a.out
#SBATCH -e ./logs/ccd_lmdb_%a.err
#SBATCH --array=0

CONFIG_PATH="configs/ccd_lmdb.yaml"
MAP_SIZE=1000000000000 # ~1TB
N_SHARDS=1
export PYTHONPATH="."

python -u scripts/build_lmdb.py build \
  ${CONFIG_PATH} \
  --map-size "${MAP_SIZE}" \
  --shard-idx "${SLURM_ARRAY_TASK_ID}" \
  -n "${N_SHARDS}" \
