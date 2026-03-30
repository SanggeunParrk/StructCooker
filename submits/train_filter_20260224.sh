#!/bin/bash
#SBATCH -J train_filter2
#SBATCH --ntasks-per-node=1   
#SBATCH -c 112     
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/train_filter_20260224.out
#SBATCH -e ./logs/train_filter_20260224.err

CONFIG_PATH="configs/train_filter_20260224.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"