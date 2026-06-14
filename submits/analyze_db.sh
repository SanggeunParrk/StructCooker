#!/bin/bash
#SBATCH -J analyze_DB
#SBATCH --ntasks-per-node=1
#SBATCH -c 68
#SBATCH --mem=270g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/analyze_DB.out
#SBATCH -e ./logs/analyze_DB.err

export PYTHONPATH="."

python -m datacooker.cli.workflow extract-lmdb "configs/analyze_train_20260224.yaml"
python -m datacooker.cli.workflow extract-lmdb "configs/analyze_train.yaml"
python -m datacooker.cli.workflow extract-lmdb "configs/analyze_valid1.yaml"
python -m datacooker.cli.workflow extract-lmdb "configs/analyze_valid2.yaml"