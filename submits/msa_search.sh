#!/bin/bash
#SBATCH -J msa_search
<<<<<<< Updated upstream
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=80g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/msa_search_%j.out
#SBATCH -e ./logs/msa_search_%j.err

set -euo pipefail
=======
#SBATCH --nodes=9
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/msa_search.out
#SBATCH -e ./logs/msa_search.err
>>>>>>> Stashed changes

export PYTHONPATH="."
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

srun --ntasks="${SLURM_NTASKS}" --ntasks-per-node=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
  python -m scripts.postprocess data_transform_parallel "configs/msa_search.yaml"
