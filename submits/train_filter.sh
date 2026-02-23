#!/bin/bash
#SBATCH -J train_filter
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=800g
#SBATCH --qos cssb_h100
#SBATCH -p h100
#SBATCH -w node01
#SBATCH -o ./logs/train_filter.out
#SBATCH -e ./logs/train_filter.err

CONFIG_PATH="configs/train_filter.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"