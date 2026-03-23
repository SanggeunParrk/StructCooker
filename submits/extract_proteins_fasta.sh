#!/bin/bash

input="/data/psk6950/BioMolDB_20260224/fasta/merged.fasta"
output="/data/psk6950/BioMolDB_20260224/fasta/polypeptide_L.fasta"

awk -F'|' '
BEGIN { keep=0 }
/^>/ {
    if ($2 ~ /polypeptide\(L\)/) {
        keep=1
        print
    } else {
        keep=0
    }
    next
}
{
    if (keep) print
}
' "$input" > "$output"
