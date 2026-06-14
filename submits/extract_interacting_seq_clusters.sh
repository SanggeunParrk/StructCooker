#!/bin/bash
#SBATCH -J extract_interacting_seq_clusters
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=64g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/extract_interacting_seq_clusters.out
#SBATCH -e ./logs/extract_interacting_seq_clusters.err

CONFIG_PATH="configs/extract_interacting_seq_clusters.yaml"
export PYTHONPATH="."

python -m datacooker.cli.workflow run "${CONFIG_PATH}"
