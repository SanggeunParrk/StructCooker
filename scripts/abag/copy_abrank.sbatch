#!/bin/bash
#SBATCH -J abrank_cp
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH -t 01:00:00
#SBATCH -o /home/psk6950/StructCooker/logs/staging/abrank_cp_%j.out
# AbRank: copy psh's AbRank_dataset_updated_v2.csv.bak_preCluster (the table before psh added the Ab_CDR90_cluster and
# Ag40_cluster columns) from cssb-master. The sha256 on both ends must be equal; a mismatch exits non-zero.
set -euo pipefail
SRC=/home/psh/data/affinity/abrank/AbRank_dataset_updated_v2.csv.bak_preCluster
R=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}/materials/raw/abrank
mkdir -p $R
cd $R
want=$(ssh -o BatchMode=yes cssb "sha256sum $SRC" | cut -d' ' -f1)
echo "source sha256 $want"
ssh -o BatchMode=yes cssb "stat -c '%s %y' $SRC"
rsync -t --partial -e "ssh -o BatchMode=yes" cssb:$SRC AbRank_dataset_updated_v2.csv.part
got=$(sha256sum AbRank_dataset_updated_v2.csv.part | cut -d' ' -f1)
[ "$got" = "$want" ] || { echo "SHA256 MISMATCH $got != $want"; exit 1; }
mv AbRank_dataset_updated_v2.csv.part AbRank_dataset_updated_v2.csv
chmod a-w AbRank_dataset_updated_v2.csv
printf 'file\tbytes\tsha256\tsource\n' > SOURCE.tsv
printf 'AbRank_dataset_updated_v2.csv\t%s\t%s\tcssb-master:%s\n' "$(stat -c %s AbRank_dataset_updated_v2.csv)" "$got" "$SRC" >> SOURCE.tsv
echo done
