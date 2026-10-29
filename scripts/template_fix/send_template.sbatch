#!/bin/bash
#SBATCH -J send_tmpl
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 4
#SBATCH --mem=8g
#SBATCH -t 12:00:00
#SBATCH -o logs/template_fix/send_template_%j.out
# template fix (2026-10-29): chain/template/template.lmdb (+17,356 chains) and lmdb/foldbench/template to H100 and cssb.
# The 154 GB data.mdb changed in place, so it is sent as a delta (--inplace --no-whole-file); rerun the script if it stops midway.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
R=/data/shared/cssb_data/BioMol_clean
KEY=$HOME/.ssh/id_ed25519
send() {
  T=$1
  if [ $T = h100 ]; then RM=psk6950@147.47.135.34; P=7777; D=data/BioMolDB_20261014; else RM=psk6950@147.47.218.30; P=7994; D=/public_data02/BioMolDB_20261014; fi
  SSH="ssh -o BatchMode=yes -p $P -i $KEY"
  cd $R
  rsync -aR --inplace --no-whole-file --exclude=lock.mdb --info=stats1 -e "$SSH" \
    lmdb/chain/template/template.lmdb lmdb/chain/template/template.lmdb.index.tsv lmdb/chain/template/template.lmdb.meta.json \
    lmdb/foldbench/template lmdb/foldbench/foldbench_build.json $RM:$D/ || { echo "RSYNC FAILED $T"; return; }
  FILES="lmdb/chain/template/template.lmdb lmdb/chain/template/template.lmdb.index.tsv lmdb/chain/template/template.lmdb.meta.json lmdb/foldbench/template lmdb/foldbench/foldbench_build.json"
  find $FILES -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort > /tmp/st_l.$T.$$
  $SSH $RM "cd $D && find $FILES -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort" > /tmp/st_r.$T.$$
  if diff /tmp/st_l.$T.$$ /tmp/st_r.$T.$$ >/dev/null; then echo "$T SIZES OK: $(wc -l < /tmp/st_l.$T.$$) files"; else echo "$T MISMATCH"; diff /tmp/st_l.$T.$$ /tmp/st_r.$T.$$ | head; fi
  for f in lmdb/chain/template/template.lmdb.index.tsv lmdb/chain/template/template.lmdb.meta.json lmdb/foldbench/foldbench_build.json lmdb/foldbench/template/template.lmdb/data.mdb; do
    a=$(md5sum $R/$f | cut -c1-12); b=$($SSH $RM "md5sum $D/$f" | cut -c1-12); [ "$a" = "$b" ] && echo "$T md5 ok $f" || echo "$T MD5 MISMATCH $f"; done
}
send h100 &
send cssb &
wait
