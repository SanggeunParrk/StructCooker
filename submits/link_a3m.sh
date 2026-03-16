#!/bin/bash
#SBATCH -J link_a3m
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/link_a3m.out
#SBATCH -e ./logs/link_a3m.err


SRC_DIR="/data/psk6950/msa"
DEST_DIR="/data/psk6950/template_db/pdb100_20260224/a3m"

find "$SRC_DIR" -type f -name "P*.a3m" | while read -r file; do
    ln -s "$file" "$DEST_DIR/$(basename "$file")"
done
