#!/bin/bash
#SBATCH -J attach_seq_metadata
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=800g
#SBATCH --qos cssb_h100
#SBATCH -p h100
#SBATCH -w node01
#SBATCH -o ./logs/attach_seq_metadata.out
#SBATCH -e ./logs/attach_seq_metadata.err

export PYTHONPATH="."

python -m scripts.build_lmdb rebuild "configs/attach_seq_metadata_valid_1.yaml"
python -m scripts.build_lmdb rebuild "configs/attach_seq_metadata_train.yaml"