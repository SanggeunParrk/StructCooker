#!/bin/bash
#SBATCH -J send_cssb3
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -c 8
#SBATCH --mem=16g
#SBATCH --exclude=node01,node04,node07
#SBATCH -t 2-00:00:00
#SBATCH -o /home/psk6950/StructCooker/logs/transfer/send_cssb3_%j.out
# 2026-11-01: bring cssb3 (210.109.54.14:/public_data/BioMolDB_20261014, last sent 10-15) to the cssb set.
#   stale -> refreshed in place : chain template, cif_pdb_attached (all PDB, train_20210930), valid_1/2, docs
#   new                         : interface_pdb, foldbench, ccd_inference, metadata (ab_seq_id, seq_cluster40_PDB, train_item, resources)
#   removed                     : metadata/seq_cluster30_PDB.tsv (30% clusters are deleted everywhere)
# One rsync stream per unit; then every file of every unit is compared: size, full md5 under 200 MB, else 8 random 64 MiB regions.
set -uo pipefail
R=/data/shared/cssb_data/BioMol_clean
KEY=$HOME/.ssh/id_ed25519
RM=psk6950@210.109.54.14; P=7994; D=/public_data/BioMolDB_20261014
SSH_E="ssh -o BatchMode=yes -p $P -i $KEY -c aes128-gcm@openssh.com"
SSH="ssh -n -o BatchMode=yes -p $P -i $KEY -c aes128-gcm@openssh.com"
UNITS="lmdb/chain/template/template.lmdb lmdb/chain/template/template.lmdb.index.tsv lmdb/chain/template/template.lmdb.meta.json
lmdb/pdb/cif_attached/cif_pdb_attached.lmdb lmdb/pdb/cif_attached/cif_pdb_attached.lmdb.index.tsv lmdb/pdb/cif_attached/cif_pdb_attached.lmdb.meta.json
lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb.index.tsv lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb.meta.json
lmdb/pdb/valid/cif_attached_valid_1.lmdb lmdb/pdb/valid/cif_attached_valid_1.lmdb.index.tsv lmdb/pdb/valid/cif_attached_valid_1.lmdb.meta.json
lmdb/pdb/valid/cif_attached_valid_2.lmdb lmdb/pdb/valid/cif_attached_valid_2.lmdb.index.tsv lmdb/pdb/valid/cif_attached_valid_2.lmdb.meta.json
lmdb/interface/pdb
lmdb/foldbench
lmdb/ccd/ccd_inference.lmdb lmdb/ccd/ccd_inference.lmdb.index.tsv lmdb/ccd/ccd_inference.lmdb.meta.json
metadata/ab_seq_id.tsv metadata/seq_cluster40_PDB.tsv metadata/train_item_20210930.tsv metadata/resources.lmdb
docs"
# one stream per group (a DB with its index and meta travels together)
GROUPS_=(
 "lmdb/chain/template/template.lmdb lmdb/chain/template/template.lmdb.index.tsv lmdb/chain/template/template.lmdb.meta.json"
 "lmdb/pdb/cif_attached/cif_pdb_attached.lmdb lmdb/pdb/cif_attached/cif_pdb_attached.lmdb.index.tsv lmdb/pdb/cif_attached/cif_pdb_attached.lmdb.meta.json"
 "lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb.index.tsv lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb.meta.json"
 "lmdb/pdb/valid/cif_attached_valid_1.lmdb lmdb/pdb/valid/cif_attached_valid_1.lmdb.index.tsv lmdb/pdb/valid/cif_attached_valid_1.lmdb.meta.json lmdb/pdb/valid/cif_attached_valid_2.lmdb lmdb/pdb/valid/cif_attached_valid_2.lmdb.index.tsv lmdb/pdb/valid/cif_attached_valid_2.lmdb.meta.json"
 "lmdb/interface/pdb"
 "lmdb/foldbench lmdb/ccd/ccd_inference.lmdb lmdb/ccd/ccd_inference.lmdb.index.tsv lmdb/ccd/ccd_inference.lmdb.meta.json"
 "metadata/ab_seq_id.tsv metadata/seq_cluster40_PDB.tsv metadata/train_item_20210930.tsv metadata/resources.lmdb docs"
)
cd $R
$SSH $RM "mkdir -p $D"
send_group() {
  local g="$1" try
  for try in 1 2 3 4; do
    rsync -aR --inplace --no-whole-file --partial --exclude=lock.mdb --exclude='.*build*' -e "$SSH_E" $g $RM:$D/ && return 0
    echo "[send] try $try failed: ${g:0:80}"; sleep 60
  done
  return 1
}
pids=()
for g in "${GROUPS_[@]}"; do send_group "$g" & pids+=($!); done
fail=0; for p in "${pids[@]}"; do wait $p || fail=1; done
echo "== sent (failures: $fail); removing the 30% clustering file"
$SSH $RM "rm -f $D/metadata/seq_cluster30_PDB.tsv; ls $D/metadata"

verify_file() {  # relative path
  local f=$1 size a b blocks off i
  size=$(stat -c %s $R/$f)
  [ "$($SSH $RM "stat -c %s $D/$f" 2>/dev/null)" = "$size" ] || { echo "SIZE MISMATCH $f"; return 1; }
  if [ $size -lt 209715200 ]; then
    [ "$(md5sum $R/$f | cut -c1-32)" = "$($SSH $RM "md5sum $D/$f" | cut -c1-32)" ] || { echo "MD5 MISMATCH $f"; return 1; }
  else
    blocks=$(( size / 67108864 )); [ $blocks -lt 1 ] && blocks=1
    for i in 0 1 2 3 4 5 6 7; do
      off=$(( (i * 2654435761 + size) % blocks ))
      a=$(dd if=$R/$f bs=64M skip=$off count=1 2>/dev/null | md5sum | cut -c1-32)
      b=$($SSH $RM "dd if=$D/$f bs=64M skip=$off count=1 2>/dev/null | md5sum" | cut -c1-32)
      [ "$a" = "$b" ] || { echo "REGION MISMATCH $f block $off"; return 1; }
    done
  fi
}
export -f verify_file; export R D RM SSH
echo "== verify"
find $UNITS -type f ! -name lock.mdb ! -path '*/.*build*' ! -name '*.pre_template_fix' > /tmp/cssb3_units.$$ 2>/dev/null
wc -l < /tmp/cssb3_units.$$
bad=$(xargs -a /tmp/cssb3_units.$$ -P 6 -I{} bash -c 'verify_file {}' | tee /dev/stderr | wc -l)
# nothing extra left in the units' directories
echo "verify problems: $bad"
echo "ALL DONE (send failures: $fail, verify problems: $bad)"
