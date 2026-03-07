#!/bin/bash
#SBATCH -J extract_fasta
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=400g
#SBATCH --qos cssb_h100
#SBATCH -p h100
#SBATCH -w node02
#SBATCH -o ./logs/extract_fasta.out
#SBATCH -e ./logs/extract_fasta.err

export PYTHONPATH="."

python -m scripts.postprocess data_transform "configs/extract_ccd_cif.yaml"