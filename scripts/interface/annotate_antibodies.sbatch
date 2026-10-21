#!/bin/bash
#SBATCH -J if_antibody
#SBATCH -c 96
#SBATCH --mem=128g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 4:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/interface/antibody_%j.out
# ANARCI over every unique PDB protein sequence -> $BIOMOL_ROOT/materials/intermediate/interface/antibody_cache.tsv
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
ROOT=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}
mkdir -p "$ROOT/materials/intermediate/interface"
pixi run python -I scripts/interface/annotate_antibodies.py "$ROOT/materials/intermediate/interface/antibody_cache.tsv" 90 2>&1 | grep -v "Limiting hmmer"
