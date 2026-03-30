#!/bin/bash
#SBATCH -J merge__db
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/merge__db.out
#SBATCH -e ./logs/merge_db.err

export PYTHONPATH="."
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"


# INPUT_PATH="/data/psk6950/BioMolDB_20260224/template_shard*.lmdb"
# OUTPUT_PATH="/data/psk6950/BioMolDB_20260224/template.lmdb"

INPUT_PATH="/data/psk6950/BioMolDB_20260224/cif_shard*.lmdb"
OUTPUT_PATH="/data/psk6950/BioMolDB_20260224/cif.lmdb"

srun --ntasks="${SLURM_NTASKS}" \
     --ntasks-per-node=1 \
     --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
     python scripts/build_lmdb.py merge "${INPUT_PATH}" -o "${OUTPUT_PATH}"
