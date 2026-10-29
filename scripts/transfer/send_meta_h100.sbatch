#!/bin/bash
#SBATCH -J send_meta
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 4
#SBATCH --mem=8g
#SBATCH -o /home/psk6950/StructCooker/logs/transfer/send_meta_%j.out
set -uo pipefail
M=/data/shared/cssb_data/BioMol_clean/metadata
SSH="ssh -o BatchMode=yes -p 7777 -i $HOME/.ssh/id_ed25519 -c aes128-gcm@openssh.com"
cd $M
FILES="seq_cluster40_AFM.tsv seq_cluster40_TDM.tsv interacting_seq_clusters.tsv interacting_seq_ids.tsv seq_id_map.tsv"
rsync -a -e "$SSH" $FILES psk6950@147.47.135.34:data/BioMolDB_20261014/metadata/
for f in $FILES seq_cluster40_PDB.tsv seq_cluster40_OFD.tsv cif_metadata.tsv train_item_20210930.tsv; do
  l=$(stat -c %s $f); r=$($SSH psk6950@147.47.135.34 "stat -c %s data/BioMolDB_20261014/metadata/$f")
  [ "$l" = "$r" ] && echo "ok $f $l" || echo "MISMATCH $f $l $r"
done
