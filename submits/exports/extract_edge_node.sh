#!/bin/bash
#SBATCH -J extract_edge_node
#SBATCH --ntasks-per-node=1
#SBATCH -c 78
#SBATCH --mem=400g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/extract_edge_node.out
#SBATCH -e ./logs/extract_edge_node.err

export PYTHONPATH="."

python -m datacooker.cli.workflow extract-lmdb "configs/exports/extract_edge_node_disordered.yaml"
python -m datacooker.cli.workflow extract-lmdb "configs/exports/extract_edge_node_train.yaml"
python -m datacooker.cli.workflow extract-lmdb "configs/exports/extract_edge_node_valid1.yaml"
python -m datacooker.cli.workflow extract-lmdb "configs/exports/extract_edge_node_valid2.yaml"
python -m datacooker.cli.workflow extract-lmdb "configs/exports/extract_edge_node_train_20260224.yaml"
