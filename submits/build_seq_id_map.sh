#!/bin/bash
#SBATCH -J build_seq_id_map
#SBATCH --ntasks-per-node=1
#SBATCH -c 40
#SBATCH --mem=200g
#SBATCH -p gpu
#SBATCH -w gpu04
#SBATCH -o ./logs/build_seq_id_map.out
#SBATCH -e ./logs/build_seq_id_map.err

CONFIG_PATH="configs/build_seq_id_map.yaml"
export PYTHONPATH="."

python -m datacooker.cli.workflow run "${CONFIG_PATH}"
