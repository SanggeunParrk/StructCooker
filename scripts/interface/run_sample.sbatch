#!/bin/bash
#SBATCH -J if_sample
#SBATCH -c 2
#SBATCH --mem=32g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 0:30:00
#SBATCH -o logs/interface/sample_%j.out
# usage: sbatch run_sample.sbatch <key> [<key> ...]
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
D=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}
pixi run python -I scripts/interface/run_sample.py $D/lmdb/pdb/cif/cif_pdb.lmdb - "$@" 2>&1
