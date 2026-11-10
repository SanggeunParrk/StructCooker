#!/bin/bash
#SBATCH -J unit_census
#SBATCH -c 2
#SBATCH --mem=32g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 4:00:00
#SBATCH -o logs/unit/census_%j.out
# usage: sbatch monomer_census.sbatch <out.tsv> (--random N SEED | --ids ID ...) [--max-mb MB]
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
D=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}
P=$D/lmdb/pdb/cif
out=$1; shift
pixi run python -I scripts/unit/monomer_census.py $P/cif_pdb.lmdb $P/cif_pdb.lmdb.index.tsv "$out" "$@" 2>&1
