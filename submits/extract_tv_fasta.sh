#!/bin/bash
#SBATCH -J extract_tv_fasta
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=400g
#SBATCH --qos cssb_h100
#SBATCH -p h100
#SBATCH -w node02
#SBATCH -o ./logs/extract_tv_fasta.out
#SBATCH -e ./logs/extract_tv_fasta.err

export PYTHONPATH="."

python -m scripts.postprocess db_extract "configs/extract_fasta_train.yaml"
python -m scripts.postprocess db_extract "configs/extract_fasta_valid_1.yaml"