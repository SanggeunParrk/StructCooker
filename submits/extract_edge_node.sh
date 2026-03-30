#!/bin/bash
#SBATCH -J extract_edge_node
#SBATCH --ntasks-per-node=1
#SBATCH -c 80
#SBATCH --mem=400g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/extract_edge_node.out
#SBATCH -e ./logs/extract_edge_node.err

export PYTHONPATH="."

python -m scripts.postprocess db_extract "configs/extract_edge_node_train.yaml"
python -m scripts.postprocess db_extract "configs/extract_edge_node_valid1.yaml"
python -m scripts.postprocess db_extract "configs/extract_edge_node_valid2.yaml"