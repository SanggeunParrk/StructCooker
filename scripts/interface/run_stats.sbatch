#!/bin/bash
#SBATCH -J if_stats
#SBATCH -c 2
#SBATCH --mem=32g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 4:00:00
#SBATCH -o logs/interface/stats_%j.out
# usage: sbatch run_stats.sbatch <N> <seed> <out_dir>
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
D=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}
P=$D/lmdb/pdb/cif
pixi run python -I scripts/interface/run_stats.py $P/cif_pdb.lmdb $P/cif_pdb.lmdb.index.tsv "$@" 2>&1
