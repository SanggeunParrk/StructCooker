#!/bin/bash
#SBATCH -J if_full
#SBATCH -c 112
#SBATCH --mem=480g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 12:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/interface/full_%j.out
# Whole PDB, three size tiers in one allocation (resumable: rerun the same command to continue).
# Peak memory is ~16x the compressed record (5j7v, 1.3 GB -> 21 GB), so the tiers trade workers for size.
# Submit from the repository root with logs/interface/ present: `structcooker interface build` does both.
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
ROOT=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}
OUT=${1:-$ROOT/lmdb/interface/pdb}
AB=$ROOT/materials/intermediate/interface/antibody_cache.tsv
run() { pixi run python -I scripts/interface/build_db.py --out "$OUT" --antibody-cache "$AB" "$@" 2>&1 | grep -v -e "Limiting hmmer" -e WARN; }
run --workers 104 --min-mb 0   --max-mb 30   --chunk 8
run --workers 24  --min-mb 30  --max-mb 300  --chunk 1
run --workers 3   --min-mb 300 --max-mb 1e6  --chunk 1
