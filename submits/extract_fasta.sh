#!/bin/bash
#SBATCH -J extract_fasta
#SBATCH --ntasks-per-node=1
#SBATCH -c 40
#SBATCH --mem=200g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/extract_fasta.out
#SBATCH -e ./logs/extract_fasta.err

export PYTHONPATH="."

python -m scripts.postprocess db_extract "configs/extract_fasta_whole.yaml"