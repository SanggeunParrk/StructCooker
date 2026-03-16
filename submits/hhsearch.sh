#!/bin/bash
#SBATCH -J hhsearch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 112
#SBATCH --mem=491g
#SBATCH -p cpu-standard
#SBATCH --qos=cpu-standard-q
#SBATCH -o ./logs/hhsearch.out
#SBATCH -e ./logs/hhsearch.err

hhsearch -b 50 -B 500 -z 50 -Z 500 -mact 0.05 -cpu 4 -maxmem 16 -aliw 100000 -e 100 -p 5.0 -d "/data/psk6950/template_db/pdb100_20260224/" \
	-i "/data/psk6950/msa/P000/P0003600/P0003600.a3m" -o "/data/psk6950/msa/P000/P0003600/P0003600.hhr" -atab "/data/psk6950/msa/P000/P0003600/P0003600.atam" -v 0
