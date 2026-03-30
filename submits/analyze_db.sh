#!/bin/bash
#SBATCH -J analyze_DB
#SBATCH --ntasks-per-node=1
#SBATCH -c 80
#SBATCH --mem=400g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/analyze_DB.out
#SBATCH -e ./logs/analyze_DB.err

export PYTHONPATH="."

python -m scripts.postprocess db_extract "configs/analyze_train.yaml"
python -m scripts.postprocess db_extract "configs/analyze_valid1.yaml"
python -m scripts.postprocess db_extract "configs/analyze_valid2.yaml"