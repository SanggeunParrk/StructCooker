#!/bin/bash
#SBATCH -J disordered_data_filter
#SBATCH --ntasks-per-node=1
#SBATCH -c 72
#SBATCH --mem=752g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/disordered_data_filter.out
#SBATCH -e ./logs/disordered_data_filter.err

CONFIG_PATH="configs/disordered_data_filter.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"