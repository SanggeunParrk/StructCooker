#!/bin/bash
#SBATCH -J attach_seq_metadata
#SBATCH --ntasks-per-node=1
#SBATCH -c 48
#SBATCH --mem=376g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/attach_seq_metadata.out
#SBATCH -e ./logs/attach_seq_metadata.err

export PYTHONPATH="."

python -m datacooker.cli.lmdb rebuild "configs/metadata/attach_seq_metadata_disordered.yaml"
python -m datacooker.cli.lmdb rebuild "configs/metadata/attach_seq_metadata_valid_1.yaml"
python -m datacooker.cli.lmdb rebuild "configs/metadata/attach_seq_metadata_train.yaml"
python -m datacooker.cli.lmdb rebuild "configs/metadata/attach_seq_metadata_train_20260224.yaml"
