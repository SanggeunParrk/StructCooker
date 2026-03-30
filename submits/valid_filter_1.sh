#!/bin/bash
#SBATCH -J valid_filter_1
#SBATCH --ntasks-per-node=1         
#SBATCH -c 112     
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/valid_filter_1.out
#SBATCH -e ./logs/valid_filter_1.err

CONFIG_PATH="configs/valid_filter_1.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"