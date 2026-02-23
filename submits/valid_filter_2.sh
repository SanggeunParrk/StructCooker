#!/bin/bash
#SBATCH -J valid_filter_2
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=800g
#SBATCH --qos cssb_h100
#SBATCH -p h100
#SBATCH -w node02
#SBATCH -o ./logs/valid_filter_2.out
#SBATCH -e ./logs/valid_filter_2.err

CONFIG_PATH="configs/valid_filter_2.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"