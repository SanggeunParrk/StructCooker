#!/bin/bash
#SBATCH -J hhmake
#SBATCH --nodes=9
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/hhmake.out
#SBATCH -e ./logs/hhmake.err

export PYTHONPATH="."
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

srun --ntasks="${SLURM_NTASKS}" --ntasks-per-node=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
  python -m scripts.postprocess data_transform_parallel "configs/hhmake.yaml"
