#!/bin/bash
#SBATCH -J valid_filter_2
#SBATCH --ntasks-per-node=1     
#SBATCH -c 32
#SBATCH --mem=128g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/valid_filter_2.out
#SBATCH -e ./logs/valid_filter_2.err

CONFIG_PATH="configs/filters/valid_filter_2.yaml"
export PYTHONPATH="."

python -m datacooker.cli.lmdb rebuild "${CONFIG_PATH}"