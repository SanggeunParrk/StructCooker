#!/bin/bash
#SBATCH -J if_cssb_transfer
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -c 4
#SBATCH --mem=8g
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/interface/transfer_%j.out
# lmdb/interface/pdb (the interface-unit DB, ~164 GB) -> cssb:/public_data02/BioMolDB_20261014, relative path
# kept (that tree mirrors BioMol_clean). One rsync stream, shard by shard with 2 s between connections, then a
# same-size listing and md5 of the small files on both ends. Rerun to resume (rsync --partial).
# Override with REMOTE=... DEST=... SRC=... in the environment.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
SRC=${SRC:-/data/shared/cssb_data/BioMol_clean}
DEST=${DEST:-/public_data02/BioMolDB_20261014}
REMOTE=${REMOTE:-psk6950@147.47.218.30}
SSH="ssh -o BatchMode=yes -p 7994 -i $HOME/.ssh/id_ed25519 -c aes128-gcm@openssh.com"
P=lmdb/interface/pdb
cd "$SRC"
$SSH "$REMOTE" "test -d $DEST && mkdir -p $DEST/$P"
echo "== small files $(date +%T)"
rsync -aR --partial --exclude=lock.mdb --exclude='shards/' -e "$SSH" ./$P ./$P/ "$REMOTE:$DEST/" 2>&1 | tail -3 || true
rsync -aR --partial --exclude=lock.mdb --exclude='shards/' -e "$SSH" ./$P/ "$REMOTE:$DEST/"
ls -d $P/interface_pdb.lmdb/shards/shard_*.lmdb | while read -r u; do echo "$(du -s --apparent-size -B1 "$u" | cut -f1) $u"; done \
  | sort -rn | cut -d' ' -f2 > /tmp/if_units.$SLURM_JOB_ID
n=$(wc -l < /tmp/if_units.$SLURM_JOB_ID); i=0
while read -r u; do
  i=$((i+1)); sleep 2
  rsync -aR --partial --exclude=lock.mdb -e "$SSH" "./$u" "$REMOTE:$DEST/"
  echo "sent $i/$n $u $(date +%T)"
done < /tmp/if_units.$SLURM_JOB_ID
echo "== verify: every file, same size on both ends"
find $P -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort > /tmp/if_local.$SLURM_JOB_ID
$SSH "$REMOTE" "cd $DEST && find $P -type f ! -name lock.mdb -printf '%p %s\n' | LC_ALL=C sort" > /tmp/if_remote.$SLURM_JOB_ID
diff /tmp/if_local.$SLURM_JOB_ID /tmp/if_remote.$SLURM_JOB_ID && echo "SIZES OK: $(wc -l < /tmp/if_local.$SLURM_JOB_ID) files, $(awk '{s+=$2}END{printf "%.1f GB", s/1e9}' /tmp/if_local.$SLURM_JOB_ID)"
for f in manifest.tsv manifest_full.tsv split.tsv report.md split_summary.md done.tsv interface_pdb.lmdb.meta.json interface_pdb.lmdb/shards.json; do
  l=$(md5sum $P/$f | cut -d' ' -f1); r=$($SSH "$REMOTE" "md5sum $DEST/$P/$f" | cut -d' ' -f1)
  [ "$l" = "$r" ] && echo "md5 ok $f" || { echo "MD5 MISMATCH $f"; exit 1; }
done
echo "VERIFY OK"
$SSH "$REMOTE" "du -sh $DEST/$P; df -h /public_data02 | tail -1"
