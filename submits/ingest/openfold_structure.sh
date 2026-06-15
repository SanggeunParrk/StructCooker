#!/bin/bash
#SBATCH -J of3_structure_DB
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/openfold_structure_%a.out
#SBATCH -e ./logs/openfold_structure_%a.err
#SBATCH --array=0

# Usage: sbatch submits/ingest/openfold_structure.sh <long|short|rna>
SET="${1:-long}"
CONFIG_PATH="configs/ingest/openfold_structure_${SET}.yaml"
MAP_SIZE=2000000000000 # ~2TB
N_SHARDS=1
export PYTHONPATH="."

pixi run python -u -m datacooker.cli.lmdb build \
  "${CONFIG_PATH}" \
  --map-size "${MAP_SIZE}" \
  --shard-idx "${SLURM_ARRAY_TASK_ID}" \
  --n-shards "${N_SHARDS}"
