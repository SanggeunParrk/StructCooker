#!/bin/bash
#SBATCH -J extract_interacting_seq_ids
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=324g
#SBATCH --qos cssb_h100
#SBATCH -p h100
#SBATCH -w node02
#SBATCH -o ./logs/extract_interacting_seq_ids.out
#SBATCH -e ./logs/extract_interacting_seq_ids.err

CONFIG_PATH="configs/extract_interacting_seq_ids.yaml"
export PYTHONPATH="."

python -m scripts.postprocess db_extract "${CONFIG_PATH}"
