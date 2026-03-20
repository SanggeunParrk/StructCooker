#!/bin/bash
#SBATCH -J remove_lower_from_a3m
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/remove_lower_from_a3m.out
#SBATCH -e ./logs/remove_lower_from_a3m.err

export PYTHONPATH="."
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

srun --ntasks="${SLURM_NTASKS}" --ntasks-per-node=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
  python -m scripts.postprocess data_transform_parallel "configs/remove_lower_from_a3m.yaml"
