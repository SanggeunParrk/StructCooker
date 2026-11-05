#!/bin/bash
#SBATCH -J asd_unpack
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 2
#SBATCH --mem=8g
#SBATCH -t 02:00:00
#SBATCH -o /home/psk6950/StructCooker/logs/staging/asd_unpack_%j.out
# Test the four Drive zips of ASD v2.0 (download/), then extract them into download/extract/ (merged asd_2.0/).
set -uo pipefail
R=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}/materials/raw/asd
cd $R/download
rc=0
for z in asd_2.0-*.zip; do unzip -tq $z || { echo "FAIL $z"; rc=1; }; done
[ $rc = 0 ] || exit 1
mkdir -p extract
for z in asd_2.0-*.zip; do unzip -qn $z -d extract || { echo "EXTRACT FAIL $z"; rc=1; }; done
find extract -type f | wc -l
du -sh extract
exit $rc
