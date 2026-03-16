#!/bin/bash
#SBATCH -J ffindex
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/ffindex.out
#SBATCH -e ./logs/ffindex.err


BASE_DIR="/data/psk6950/template_db/pdb100_20260224"
DB_PREFIX="${BASE_DIR}/pdb100_20260224"

HHM_SRC_DIR="/data/psk6950/hhm"
A3M_SRC_DIR="${BASE_DIR}/a3m"

HHM_FFDATA="${DB_PREFIX}_hhm.ffdata"
HHM_FFINDEX="${DB_PREFIX}_hhm.ffindex"

A3M_FFDATA="${DB_PREFIX}_a3m.ffdata"
A3M_FFINDEX="${DB_PREFIX}_a3m.ffindex"

CS219_FFDATA="${DB_PREFIX}_cs219.ffdata"
CS219_FFINDEX="${DB_PREFIX}_cs219.ffindex"

echo "=== Step 1: build HHM ffindex DB ==="
if [[ -f "$HHM_FFDATA" && -f "$HHM_FFINDEX" ]]; then
    echo "HHM DB already exists. Skipping."
else
    echo "Building HHM DB from: $HHM_SRC_DIR"
    ffindex_build -s "$HHM_FFDATA" "$HHM_FFINDEX" "$HHM_SRC_DIR"
fi

echo "=== Step 2: build A3M ffindex DB ==="
if [[ -f "$A3M_FFDATA" && -f "$A3M_FFINDEX" ]]; then
    echo "A3M DB already exists. Skipping."
else
    echo "Building A3M DB from: $A3M_SRC_DIR"
    ffindex_build -s "$A3M_FFDATA" "$A3M_FFINDEX" "$A3M_SRC_DIR"
fi

echo "=== Step 3: build CS219 DB ==="
if [[ -f "$CS219_FFDATA" && -f "$CS219_FFINDEX" ]]; then
    echo "CS219 DB already exists. Skipping."
else
    echo "Running cstranslate from A3M DB"
    cstranslate \
        -i "${DB_PREFIX}_a3m" \
        -o "${DB_PREFIX}_cs219" \
        -I a3m \
        -x 0.3 \
        -c 4 \
    	--ffindex 
fi

echo "=== Final files ==="
ls -lh \
    "$HHM_FFDATA" "$HHM_FFINDEX" \
    "$A3M_FFDATA" "$A3M_FFINDEX" \
    "$CS219_FFDATA" "$CS219_FFINDEX"

