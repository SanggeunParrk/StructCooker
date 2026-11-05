#!/bin/bash
#SBATCH -J asd_stage
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 2
#SBATCH --mem=8g
#SBATCH -t 02:00:00
#SBATCH -o /home/psk6950/StructCooker/logs/staging/asd_stage_%j.out
set -euo pipefail
python3 -I /home/psk6950/StructCooker/scripts/abag/stage_asd.py --root ${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}/materials/raw/asd
