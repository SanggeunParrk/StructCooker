#!/bin/bash
#SBATCH -J msa_remap
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=80g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -w node02
#SBATCH -o ./logs/msa_remap.out
#SBATCH -e ./logs/msa_rempa.err

python -u remap_msa.py
