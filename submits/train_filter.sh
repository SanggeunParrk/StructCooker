#!/bin/bash
#SBATCH -J train_filter
#SBATCH --ntasks-per-node=1
#SBATCH -c 32
#SBATCH --mem=360g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/train_filter.out
#SBATCH -e ./logs/train_filter.err

CONFIG_PATH="configs/train_filter.yaml"
export PYTHONPATH="."

python -m datacooker.cli.lmdb rebuild "${CONFIG_PATH}"
