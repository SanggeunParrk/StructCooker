#!/bin/bash
#SBATCH -J cif_pdb
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/%x_%A_%a.out
#SBATCH -e ./logs/%x_%A_%a.err
#SBATCH --array=0

# Usage: N_SHARDS=20 sbatch ... submits/ingest/cif_pdb.sh <light|medium|large|xlarge>
TIER="${1:?tier required}"
CONFIG_PATH="configs/ingest/cif_pdb_${TIER}.yaml"
MAP_SIZE=1000000000000  # ~1TB (sparse)
N_SHARDS="${N_SHARDS:-1}"
export PYTHONPATH="."

if [ "${N_SHARDS}" -gt 1 ]; then
  pixi run python -u -m datacooker.cli.lmdb build \
    "${CONFIG_PATH}" --map-size "${MAP_SIZE}" \
    --shard-idx "${SLURM_ARRAY_TASK_ID}" --n-shards "${N_SHARDS}"
else
  pixi run python -u -m datacooker.cli.lmdb build \
    "${CONFIG_PATH}" --map-size "${MAP_SIZE}"
fi
