#!/bin/bash
#SBATCH -J extract_query_fasta
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -o ./logs/extract_query_fasta_%j.out
#SBATCH -e ./logs/extract_query_fasta_%j.err

set -euo pipefail

INPUT_DIR="/data/psk6950/BioMolDB_20260224/a3m"
OUTPUT="/data/psk6950/BioMolDB_20260224/merged.fasta"
THREADS=112

LOCKFILE="${OUTPUT}.lock"
: > "$OUTPUT"
: > "$LOCKFILE"

export OUTPUT LOCKFILE

extract_one() {
    local zst="$1"
    local name
    name=$(basename "$zst" .a3m.zst)

    local seq
    seq=$(zstd -dcq "$zst" | awk '
        /^#/ { next }
        /^>/ {
            if (seen) exit
            seen = 1
            next
        }
        seen { printf "%s", $0 }
    ')

    if [[ -n "$seq" ]]; then
        (
            flock -x 9
            printf ">%s\n%s\n" "$name" "$seq"
        ) 9>>"$LOCKFILE" >> "$OUTPUT"
    fi
}
export -f extract_one

find "$INPUT_DIR" -type f -name '*.a3m.zst' -print0 \
    | xargs -0 -n 1 -P "$THREADS" bash -c 'extract_one "$0"'

rm -f "$LOCKFILE"
echo "Done: $(grep -c '^>' "$OUTPUT") sequences written to $OUTPUT"
