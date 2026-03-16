#!/bin/bash
#SBATCH -J valid_filter_1
#SBATCH --ntasks-per-node=1
#SBATCH -c 72
#SBATCH --mem=328g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/valid_filter_1.out
#SBATCH -e ./logs/valid_filter_1.err

CONFIG_PATH="configs/valid_filter_1.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"