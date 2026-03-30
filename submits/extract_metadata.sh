#!/bin/bash
#SBATCH -J extract_metadata
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/extract_metadata.out
#SBATCH -e ./logs/extract_metadata.err

CONFIG_PATH="configs/extract_metadata.yaml"
export PYTHONPATH="."

python -m scripts.postprocess db_extract "${CONFIG_PATH}"
