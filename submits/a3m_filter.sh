#!/bin/bash
#SBATCH -J a3m_filter
#SBATCH --ntasks-per-node=1   
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/a3m_filter.out
#SBATCH -e ./logs/a3m_filter.err

CONFIG_PATH="configs/a3m_filter.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"
