#!/bin/bash
#SBATCH -J merge_template_db
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/merge_template_db.out
#SBATCH -e ./logs/merge_template_db.err

export PYTHONPATH="."
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

srun --ntasks="${SLURM_NTASKS}" --ntasks-per-node=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
	python scripts/build_lmdb.py merge "/data/psk6950/BioMolDB_20260224/template_shard*.lmdb" -o "/data/psk6950/BioMolDB_20260224/template.lmdb"

