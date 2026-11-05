#!/bin/bash
#SBATCH -J poisk_dl
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH -t 01:00:00
#SBATCH -o /home/psk6950/StructCooker/logs/staging/poisk_dl_%j.out
# POISK (Zenodo record 19883049, CC BY 4.0): the four files of the record, each checked against the md5 Zenodo lists.
# Any download error or md5 mismatch exits non-zero (an earlier gdown run of ASD exited 0 with nothing downloaded).
set -euo pipefail
REC=19883049
R=${BIOMOL_ROOT:-/data/shared/cssb_data/BioMol_clean}/materials/raw/poisk/zenodo_$REC
mkdir -p $R
cd $R
curl -fsS --retry 3 --max-time 60 https://zenodo.org/api/records/$REC -o record.json
python3 -I - <<'PY' > files.tsv
import json
for f in json.load(open("record.json"))["files"]:
    print(f["key"], f["size"], f["checksum"].split(":", 1)[1], f["links"]["self"], sep="\t")
PY
cat files.tsv
while IFS=$'\t' read -r name size md5 url; do
    curl -fSL --retry 3 --max-time 1800 -o "$name.part" "$url"
    got=$(md5sum "$name.part" | cut -d' ' -f1)
    [ "$got" = "$md5" ] || { echo "MD5 MISMATCH $name: $got != $md5"; exit 1; }
    [ "$(stat -c %s "$name.part")" = "$size" ] || { echo "SIZE MISMATCH $name"; exit 1; }
    mv "$name.part" "$name"
    echo "ok $name $size $md5"
done < files.tsv
printf 'file\tbytes\tmd5\tsha256\n' > SOURCE.tsv
while IFS=$'\t' read -r name size md5 url; do
    printf '%s\t%s\t%s\t%s\n' "$name" "$size" "$md5" "$(sha256sum "$name" | cut -d' ' -f1)" >> SOURCE.tsv
done < files.tsv
chmod a-w *.csv *.zip *.json
echo done
