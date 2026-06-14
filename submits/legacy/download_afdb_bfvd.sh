#!/bin/bash
#SBATCH -J dw_af_bfvd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/dw_af_bfvd.out
#SBATCH -e ./logs/dw_af_bfvd.err

lftp -e "mirror --parallel=8 -c /pub/databases/alphafold/collaborations/bfvd/ /data/psk6950/AFDB/bfvd/; quit" https://ftp.ebi.ac.uk
