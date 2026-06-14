#!/bin/bash
#SBATCH -J rename_a3m
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/rename_a3m.out
#SBATCH -e ./logs/rename_a3m.err

mkdir -p ./logs

BASE_DIR="/data/psk6950/msa"

find "$BASE_DIR" -type f -name "t000_msa0.a3m" | while read src; do
    dir=$(dirname "$src")
    pid=$(basename "$dir")
    dst="${dir}/${pid}.a3m"

    mv "$src" "$dst"
    echo "Renamed: $src -> $dst"
done
