#!/bin/bash
#SBATCH -J train_filter
#SBATCH --ntasks-per-node=1
#SBATCH -c 64
#SBATCH --mem=600g
#SBATCH -p gpu
#SBATCH -w gpu04
#SBATCH -o ./logs/train_filter.out
#SBATCH -e ./logs/train_filter.err

CONFIG_PATH="configs/train_filter.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"