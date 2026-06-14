#!/bin/bash
#SBATCH -J train_filter2
#SBATCH --ntasks-per-node=1   
#SBATCH -c 32
#SBATCH --mem=360g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/train_filter_20260224.out
#SBATCH -e ./logs/train_filter_20260224.err

CONFIG_PATH="configs/filters/train_filter_20260224.yaml"
export PYTHONPATH="."

python -m datacooker.cli.lmdb rebuild "${CONFIG_PATH}"