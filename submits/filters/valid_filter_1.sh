#!/bin/bash
#SBATCH -J valid_filter_1
#SBATCH --ntasks-per-node=1         
#SBATCH -c 32
#SBATCH --mem=128g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/valid_filter_1.out
#SBATCH -e ./logs/valid_filter_1.err

CONFIG_PATH="configs/filters/valid_filter_1.yaml"
export PYTHONPATH="."

python -m datacooker.cli.lmdb rebuild "${CONFIG_PATH}"