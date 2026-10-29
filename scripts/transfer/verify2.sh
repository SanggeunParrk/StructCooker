#!/bin/bash
# verify2.sh h100|cssb [group1only] -- every file of files_<target>.tsv present on the remote with the same size
T=$1
if [ "${T#h100}" != "$T" ]; then REMOTE=psk6950@147.47.135.34; PORT=7777; DEST=data/BioMolDB_20261014; else REMOTE=psk6950@147.47.218.30; PORT=7994; DEST=/public_data02/BioMolDB_20261014; fi
HERE=/home/psk6950/StructCooker/logs/transfer
LC_ALL=C sort $HERE/files_$T.tsv | tr '\t' ' ' > /tmp/v2_local.$$
cut -d' ' -f1 /tmp/v2_local.$$ | ssh -o BatchMode=yes -p $PORT -i $HOME/.ssh/id_ed25519 $REMOTE "cd $DEST && xargs -d '\n' stat -c '%n %s' 2>/dev/null | LC_ALL=C sort" > /tmp/v2_remote.$$
echo "local $(wc -l < /tmp/v2_local.$$) remote $(wc -l < /tmp/v2_remote.$$)"
comm -23 /tmp/v2_local.$$ /tmp/v2_remote.$$ | head -${N:-10} > /tmp/v2_diff.$$
if [ -s /tmp/v2_diff.$$ ]; then echo "MISSING/DIFFERENT (first ${N:-10}):"; cat /tmp/v2_diff.$$; echo "count: $(comm -23 /tmp/v2_local.$$ /tmp/v2_remote.$$ | wc -l)"; else echo "VERIFY OK"; fi
