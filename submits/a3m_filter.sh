#!/bin/bash
#SBATCH -J a3m_filter
#SBATCH --ntasks-per-node=1   
#SBATCH -c 80
#SBATCH --mem=488g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/a3m_filter_4k.out
#SBATCH -e ./logs/a3m_filter_4k.err

CONFIG_PATH="configs/a3m_filter.yaml"
export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "${CONFIG_PATH}"
