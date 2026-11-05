#!/bin/bash
#SBATCH -J antibinder_dl
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH -t 01:00:00
#SBATCH -o /home/psk6950/StructCooker/logs/staging/antibinder_dl_%j.out
# AntiBinder (Brief Bioinform 26(1) bbaf008): github.com/fdu-wangfeilab/AntiBinder, MIT licence. The original source of ASD's covid-19, hiv, biomap
# and met rows. Cloned as it is; the commit is recorded in SOURCE.tsv.
set -euo pipefail
R=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}/materials/raw/antibinder
mkdir -p $R
cd $R
[ -d repo ] && { echo "repo exists"; exit 1; }
git clone --quiet https://github.com/fdu-wangfeilab/AntiBinder repo
cd repo
commit=$(git rev-parse HEAD)
date=$(git log -1 --format=%cI)
echo "commit $commit $date"
cd ..
printf 'what\tvalue\n' > SOURCE.tsv
printf 'url\thttps://github.com/fdu-wangfeilab/AntiBinder\ncommit\t%s\ncommit_date\t%s\nfetched\t%s\n' "$commit" "$date" "$(date -I)" >> SOURCE.tsv
find repo -path repo/.git -prune -o -type f -print | xargs -d '\n' sha256sum > repo.sha256
wc -l repo.sha256
du -sh repo
echo done
