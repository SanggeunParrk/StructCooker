#!/bin/bash
#SBATCH -J build_seq_id_map
#SBATCH --ntasks-per-node=1
#SBATCH -c 32
#SBATCH --mem=64g
#SBATCH -p h100
#SBATCH --qos cssb_h100
#SBATCH -w node02
#SBATCH -o ./logs/build_seq_id_map.out
#SBATCH -e ./logs/build_seq_id_map.err

CONFIG_PATH="configs/build_seq_id_map.yaml"
export PYTHONPATH="."

python -m scripts.postprocess data_transform "${CONFIG_PATH}"
