#!/bin/bash
#SBATCH -J msa_search
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH --mem=656g
#SBATCH -p gpu
#SBATCH -w gpu05
#SBATCH -o ./logs/msa_search.out
#SBATCH -e ./logs/msa_search.err

export PYTHONPATH="."
python -m scripts.postprocess data_transform_parallel "configs/msa_search.yaml"
