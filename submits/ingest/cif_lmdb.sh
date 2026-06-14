#!/bin/bash
#SBATCH -J cif_DB_급함
#SBATCH --ntasks-per-node=1       
#SBATCH -c 96
#SBATCH --mem=720g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/cif_lmdb.out
#SBATCH -e ./logs/cif_lmdb.err

CONFIG_PATH="configs/ingest/cif_lmdb.yaml"
MAP_SIZE=1000000000000 # ~1TB
export PYTHONPATH="."

python -u -m datacooker.cli.lmdb build \
  ${CONFIG_PATH} \
  --map-size "${MAP_SIZE}" 
