#!/bin/bash
#SBATCH -J analyze_valid2
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=800g
#SBATCH --qos cssb_h100
#SBATCH -p h100
#SBATCH -w node02
#SBATCH -o ./logs/analyze_valid2.out
#SBATCH -e ./logs/analyze_valid2.err

export PYTHONPATH="."

python -m scripts.postprocess db_extract "configs/analyze_train.yaml"
python -m scripts.postprocess db_extract "configs/analyze_valid1.yaml"
python -m scripts.postprocess db_extract "configs/analyze_valid2.yaml"