#!/bin/bash
#SBATCH -J a3m_DB
#SBATCH --ntasks-per-node=1         
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/a3m_lmdb.out
#SBATCH -e ./logs/a3m_lmdb.err
#SBATCH --array=0

CONFIG_PATH="configs/a3m_lmdb.yaml"
MAP_SIZE=2000000000000 # ~2TB
N_SHARDS=1
export PYTHONPATH="."

python -u scripts/build_lmdb.py build \
  ${CONFIG_PATH} \
  --map-size "${MAP_SIZE}" \
  --shard-idx "${SLURM_ARRAY_TASK_ID}" \
  -n "${N_SHARDS}"
