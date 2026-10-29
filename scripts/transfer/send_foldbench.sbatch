#!/bin/bash
#SBATCH -J send_foldbench
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -c 4
#SBATCH --mem=8g
#SBATCH -o /home/psk6950/StructCooker/logs/transfer/send_foldbench_%j.out
# The FoldBench DB cut from the release (2026-10-29) -> H100:~/data/ and cssb:/public_data02/, same directory name.
set -uo pipefail
O=/data/psk6950/foldbench_from_clean_20261029
N=$(basename $O)
KEY=$HOME/.ssh/id_ed25519
for T in h100 cssb; do
  if [ $T = h100 ]; then R=psk6950@147.47.135.34; P=7777; D=data; else R=psk6950@147.47.218.30; P=7994; D=/public_data02; fi
  SSH="ssh -o BatchMode=yes -p $P -i $KEY"
  rsync -a --exclude=lock.mdb -e "$SSH" $O $R:$D/ || { echo "RSYNC FAILED $T"; continue; }
  (cd $O && find . -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort) > /tmp/fb_l.$$
  $SSH $R "cd $D/$N && find . -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort" > /tmp/fb_r.$$
  if diff /tmp/fb_l.$$ /tmp/fb_r.$$ >/dev/null; then echo "$T SIZES OK: $(wc -l < /tmp/fb_l.$$) files"; else echo "$T MISMATCH"; diff /tmp/fb_l.$$ /tmp/fb_r.$$ | head; fi
  for f in foldbench_build.json; do a=$(md5sum $O/$f | cut -c1-12); b=$($SSH $R "md5sum $D/$N/$f" | cut -c1-12); [ "$a" = "$b" ] && echo "$T md5 ok $f" || echo "$T MD5 MISMATCH $f"; done
done
