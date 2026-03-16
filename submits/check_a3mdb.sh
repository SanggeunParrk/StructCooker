#!/bin/bash
#SBATCH -J check_db
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/check_db.out
#SBATCH -e ./logs/check_db.err


export PYTHONPATH="."

python -m scripts.build_lmdb check_db /data/psk6950/BioMolDB_20260224/a3m.lmdb
