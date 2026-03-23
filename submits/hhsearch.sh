#!/bin/bash
#SBATCH -J hhsearch
#SBATCH --nodes=9
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/hhsearch.out
#SBATCH -e ./logs/hhsearch.err

export PYTHONPATH="."
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

srun --ntasks="${SLURM_NTASKS}" --ntasks-per-node=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
  python -m scripts.postprocess data_transform_parallel "configs/hhsearch.yaml"
