#!/bin/bash
#SBATCH -J extract_tv_fasta
#SBATCH --ntasks-per-node=1
#SBATCH -c 40
#SBATCH --mem=200g
#SBATCH -p gpu
#SBATCH -w gpu04
#SBATCH -o ./logs/extract_tv_fasta.out
#SBATCH -e ./logs/extract_tv_fasta.err

export PYTHONPATH="."

python -m datacooker.cli.workflow extract-lmdb "configs/exports/extract_fasta_train.yaml"
python -m datacooker.cli.workflow extract-lmdb "configs/exports/extract_fasta_valid_1.yaml"