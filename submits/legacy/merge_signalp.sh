#!/bin/bash
#SBATCH -J merge_signalp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/merge_signalp.out
#SBATCH -e ./logs/merge_signalp.err

SRC_ROOT="/data/psk6950/msa"
DST_ROOT="/data/psk6950/signalp"

mkdir -p "$DST_ROOT"
mkdir -p ./logs

export SRC_ROOT DST_ROOT

find "$SRC_ROOT" -type f -path "*/signalp/output.gff3" | while read -r gff; do
    if [ "$(wc -l < "$gff")" -ge 2 ]; then
        # 예: /data/psk6950/msa/P000/P0009233/signalp/output.gff3 -> PID=P0009233
        pid="$(basename "$(dirname "$(dirname "$gff")")")"
        cp "$gff" "$DST_ROOT/${pid}.gff3"
        echo "Copied: $gff -> $DST_ROOT/${pid}.gff3"
    fi
done

echo "Done."
