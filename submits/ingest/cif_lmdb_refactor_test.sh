#!/bin/bash
#SBATCH -J cif_refactor_test
#SBATCH --ntasks-per-node=1
#SBATCH -c 64
#SBATCH --mem=720g
#SBATCH -p gpu
#SBATCH -w gpu04
#SBATCH -o ./logs/cif_lmdb_refactor_test.out
#SBATCH -e ./logs/cif_lmdb_refactor_test.err

cd /home/psk6950/practice/StructCooker
export PYTHONPATH="."

CONFIG_PATH="configs/ingest/cif_lmdb_refactor_test.yaml"
MAP_SIZE=1000000000000 # ~1TB

# Run inside the pixi env (datacooker/structcooker are installed there, not in
# the base anaconda python that SLURM otherwise picks up).
/home/psk6950/.pixi/bin/pixi run python -u -m datacooker.cli.lmdb build \
  "${CONFIG_PATH}" \
  --map-size "${MAP_SIZE}"
