#!/bin/bash
#SBATCH -J filter_valdation_DB
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=800g
#SBATCH -p h100
#SBATCH -w node01
#SBATCH -o ./logs/filter_validation_DB.out
#SBATCH -e ./logs/filter_validation_DB.err

CONFIG_PATH="configs/filter_cif.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"