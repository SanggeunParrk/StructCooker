#!/bin/bash
#SBATCH -J send_upd
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 8
#SBATCH --mem=16g
#SBATCH -t 12:00:00
#SBATCH -o /home/psk6950/StructCooker/logs/transfer/send_update_%j.out
# 2026-11-01: DBs changed after the 10-28/29 refresh.
#   cssb : lmdb/pdb/cif_attached/cif_pdb_attached  (cssb held the 10-12 build; delta, the 115 GB file changed in place)
#   H100 : seqid_template_mols, OFD rna msa (full, d512, d5k), disordered msa (full, d512): changed or patched, never sent before
# Each DB is its own rsync stream; then sizes, md5 of index/meta, and md5 of 8 random 64 MiB regions of every data file.
set -uo pipefail
R=/data/shared/cssb_data/BioMol_clean
KEY=$HOME/.ssh/id_ed25519
H100_DBS="lmdb/chain/template/seqid_template_mols lmdb/distillation_rna/msa/msa_rna lmdb/distillation_rna/msa/msa_rna_d512 lmdb/distillation_rna/msa/msa_rna_d5k lmdb/distillation_disordered/msa/msa_disordered lmdb/distillation_disordered/msa/msa_disordered_d512"
CSSB_DBS="lmdb/pdb/cif_attached/cif_pdb_attached"

target() {  # name -> user@host port dest
  if [ $1 = h100 ]; then echo "psk6950@147.47.135.34 7777 data/BioMolDB_20261014"; else echo "psk6950@147.47.218.30 7994 /public_data02/BioMolDB_20261014"; fi
}
send_one() {  # target db
  read RM P D <<< "$(target $1)"; SSH="ssh -n -o BatchMode=yes -p $P -i $KEY"; db=$2
  cd $R
  for try in 1 2 3; do
    rsync -aR --inplace --no-whole-file --partial --exclude=lock.mdb -e "ssh -o BatchMode=yes -p $P -i $KEY" $db.lmdb $db.lmdb.index.tsv $db.lmdb.meta.json $RM:$D/ && break
    echo "[$1] $db rsync try $try failed"; sleep 30
  done
}
verify_one() {
  read RM P D <<< "$(target $1)"; SSH="ssh -n -o BatchMode=yes -p $P -i $KEY"; db=$2; cd $R; bad=0
  find $db.lmdb $db.lmdb.index.tsv $db.lmdb.meta.json -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort > /tmp/u_l.$1.$$
  $SSH $RM "cd $D && find $db.lmdb $db.lmdb.index.tsv $db.lmdb.meta.json -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort" > /tmp/u_r.$1.$$
  diff -q /tmp/u_l.$1.$$ /tmp/u_r.$1.$$ >/dev/null || { echo "[$1] SIZE MISMATCH $db"; diff /tmp/u_l.$1.$$ /tmp/u_r.$1.$$ | head -4; bad=1; }
  for f in $db.lmdb.index.tsv $db.lmdb.meta.json; do
    [ "$(md5sum $f | cut -c1-32)" = "$($SSH $RM "md5sum $D/$f" | cut -c1-32)" ] || { echo "[$1] MD5 MISMATCH $f"; bad=1; }
  done
  while read f size; do
    case $f in *.lmdb/*) ;; *) continue;; esac
    blocks=$(( size / 67108864 )); [ $blocks -lt 1 ] && blocks=1
    for i in 0 1 2 3 4 5 6 7; do
      off=$(( (i * 2654435761 + size) % blocks ))
      a=$(dd if=$R/$f bs=64M skip=$off count=1 2>/dev/null | md5sum | cut -c1-32)
      b=$($SSH $RM "dd if=$D/$f bs=64M skip=$off count=1 2>/dev/null | md5sum" | cut -c1-32)
      [ "$a" = "$b" ] || { echo "[$1] REGION MISMATCH $f block $off"; bad=1; }
    done
  done < /tmp/u_l.$1.$$
  [ $bad = 0 ] && echo "[$1] OK $db ($(wc -l < /tmp/u_l.$1.$$) files)"
}
pids=()
for db in $H100_DBS; do ( send_one h100 $db; verify_one h100 $db ) & pids+=($!); done
for db in $CSSB_DBS; do ( send_one cssb $db; verify_one cssb $db ) & pids+=($!); done
for p in "${pids[@]}"; do wait $p; done
echo "ALL DONE"
