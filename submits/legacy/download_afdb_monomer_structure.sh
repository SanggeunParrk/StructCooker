#!/bin/bash
#SBATCH -J dw_af_ms
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/dw_af_monomer_structure.out
#SBATCH -e ./logs/dw_af_monomer_structure.err

lftp -e "mirror --parallel=8 -c /pub/databases/alphafold/latest/ /data/psk6950/AFDB/latest; quit" https://ftp.ebi.ac.uk
