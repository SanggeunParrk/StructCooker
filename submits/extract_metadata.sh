#!/bin/bash
#SBATCH -J extract_metadata
#SBATCH --ntasks-per-node=1
#SBATCH -c 128
#SBATCH --mem=800g
#SBATCH -p h100
#SBATCH -w node02
#SBATCH -o ./logs/extract_metadata.out
#SBATCH -e ./logs/extract_metadata.err

CONFIG_PATH="configs/extract_metadata.yaml"
export PYTHONPATH="."

python -m scripts.postprocess db_extract "${CONFIG_PATH}"
