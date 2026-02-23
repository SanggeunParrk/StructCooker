#!/bin/bash
#SBATCH -J seq_cluster30
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=320g
#SBATCH --qos cssb_h100
#SBATCH -p h100
#SBATCH -w node02
#SBATCH -o ./logs/seq_cluster30.out
#SBATCH -e ./logs/seq_cluster30.err

CONFIG_PATH="configs/seq_cluster30.yaml"
export PYTHONPATH="."
export PATH="$PWD/mmseqs/bin:$PATH"

python -m scripts.postprocess data_transform "${CONFIG_PATH}"
