#!/bin/bash
#SBATCH -J sc_refresh
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -c 8
#SBATCH --mem=16g
#SBATCH --time=2-00:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o /home/psk6950/StructCooker/logs/transfer/refresh_%x_%A_%a.out
# 2026-10-28 refresh -> sbatch --export=ALL,TARGET=h100|cssb --array=0-K-1 -J name transfer2.sbatch
# TARGET h100: psk6950@147.47.135.34:7777, ~/data/BioMolDB_20261014. cssb: psk6950@147.47.218.30:7994, /public_data02/BioMolDB_20261014.
# Each rank sends units_<TARGET>_rank<R>.txt (group 1 first, then group 2) with 3 rsync streams; rerun resumes.
set -uo pipefail
HERE=/home/psk6950/StructCooker/logs/transfer
if [ "${TARGET#h100}" != "$TARGET" ]; then REMOTE=psk6950@147.47.135.34; PORT=7777; DEST=data/BioMolDB_20261014
else REMOTE=psk6950@147.47.218.30; PORT=7994; DEST=/public_data02/BioMolDB_20261014; fi
export SSH="ssh -o BatchMode=yes -p $PORT -i $HOME/.ssh/id_ed25519 -c aes128-gcm@openssh.com"
export REMOTE DEST
R=${SLURM_ARRAY_TASK_ID}
echo "$TARGET rank $R on $(hostname), $(wc -l < $HERE/units_${TARGET}_rank$R.txt) files, start $(date +%T)"
while read -r p; do $SSH -n $REMOTE "mkdir -p $DEST/$p"; done < $HERE/prefixes_$TARGET.txt
tr '\t\n' '\0\0' < $HERE/units_${TARGET}_rank$R.txt | xargs -0 -n3 -P4 $HERE/send_one.sh
echo "rank $R done $(date +%T)"
