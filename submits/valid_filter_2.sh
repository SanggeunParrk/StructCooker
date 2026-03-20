#!/bin/bash
#SBATCH -J valid_filter_2
#SBATCH --ntasks-per-node=1
#SBATCH -c 72
#SBATCH --mem=328g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/valid_filter_2.out
#SBATCH -e ./logs/valid_filter_2.err

CONFIG_PATH="configs/valid_filter_2.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"