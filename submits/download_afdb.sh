#!/bin/bash
#SBATCH -J dw_af
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/dw_af.out
#SBATCH -e ./logs/dw_af.err

lftp -e "mirror --parallel=8 -c /pub/databases/alphafold/collaborations/nvda/msas/ ./msas; quit" ftp://ftp.ebi.ac.uk