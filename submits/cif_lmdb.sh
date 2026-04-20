#!/bin/bash
#SBATCH -J cif_DB_급함
#SBATCH --ntasks-per-node=1       
#SBATCH -c 96
#SBATCH --mem=720g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/cif_lmdb.out
#SBATCH -e ./logs/cif_lmdb.err

CONFIG_PATH="configs/cif_lmdb.yaml"
MAP_SIZE=1000000000000 # ~1TB
export PYTHONPATH="."

python -u scripts/build_lmdb.py build \
  ${CONFIG_PATH} \
  --map-size "${MAP_SIZE}" 
