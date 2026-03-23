#!/bin/bash
#SBATCH -J extract_metadata
#SBATCH --ntasks-per-node=1
#SBATCH -c 88
#SBATCH --mem=628g
#SBATCH -p gpu
#SBATCH -w gpu04
#SBATCH -o ./logs/extract_metadata.out
#SBATCH -e ./logs/extract_metadata.err

CONFIG_PATH="configs/extract_metadata.yaml"
export PYTHONPATH="."

python -m scripts.postprocess db_extract "${CONFIG_PATH}"
