#!/bin/bash
#SBATCH -J fb_native
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 4
#SBATCH --mem=8g
#SBATCH -o /home/psk6950/StructCooker/logs/transfer/fb_native_%j.out
# FoldBench DB into the release layout: lmdb/foldbench/{msa,template,CCD,foldbench_build.json} next to msa_edge.lmdb.
set -uo pipefail
O=/data/psk6950/foldbench_from_clean_20261029
R=/data/shared/cssb_data/BioMol_clean/lmdb/foldbench
test -e $R/msa -o -e $R/template -o -e $R/CCD && { echo "release already has msa/template/CCD under lmdb/foldbench; stop"; exit 1; }
rsync -a --exclude=lock.mdb $O/msa $O/template $O/CCD $O/foldbench_build.json $R/ && echo "release copy done"
(cd $O && find msa template CCD foldbench_build.json -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort) > /tmp/fbn_l.$$
(cd $R && find msa template CCD foldbench_build.json -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort) > /tmp/fbn_r.$$
diff /tmp/fbn_l.$$ /tmp/fbn_r.$$ >/dev/null && echo "release SIZES OK" || echo "release MISMATCH"
KEY=$HOME/.ssh/id_ed25519
for T in h100 cssb; do
  if [ $T = h100 ]; then RM=psk6950@147.47.135.34; P=7777; D=data/BioMolDB_20261014; OLD=data/foldbench_from_clean_20261029; else RM=psk6950@147.47.218.30; P=7994; D=/public_data02/BioMolDB_20261014; OLD=/public_data02/foldbench_from_clean_20261029; fi
  SSH="ssh -o BatchMode=yes -p $P -i $KEY"
  $SSH $RM "mkdir -p $D/lmdb/foldbench && test -d $OLD && cd $OLD && for x in msa template CCD foldbench_build.json; do test -e $D/lmdb/foldbench/\$x && echo EXISTS \$x || mv \$x $D/lmdb/foldbench/; done; rmdir $OLD 2>&1 | head -1"
  rsync -a --exclude=lock.mdb -e "$SSH" $R/msa_edge.lmdb $R/msa_edge.lmdb.index.tsv $R/msa_edge.lmdb.meta.json $RM:$D/lmdb/foldbench/
  (cd $R && find . -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort) > /tmp/fbn_l.$$
  $SSH $RM "cd $D/lmdb/foldbench && find . -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort" > /tmp/fbn_r.$$
  diff /tmp/fbn_l.$$ /tmp/fbn_r.$$ >/dev/null && echo "$T SIZES OK: $(wc -l < /tmp/fbn_l.$$) files" || { echo "$T MISMATCH"; diff /tmp/fbn_l.$$ /tmp/fbn_r.$$ | head; }
done
