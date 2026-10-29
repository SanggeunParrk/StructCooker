#!/bin/bash
#SBATCH -J send_complex
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -c 8
#SBATCH --mem=16g
#SBATCH --time=1-00:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o /home/psk6950/StructCooker/logs/transfer/send_complex_%j.out
# After the PDB complex design DB is validated (318864): send it to H100 and cssb, then size-verify both ends.
set -uo pipefail
HERE=/home/psk6950/StructCooker/logs/transfer
cd /home/psk6950/StructCooker
test -f /data/psk6950/design/db_complex_pdb/manifest.tsv || { echo "no manifest.tsv"; exit 1; }
for T in h100c cssbc; do
  python3 $HERE/plan2.py $T 1
  if [ "$T" = h100c ]; then REMOTE=psk6950@147.47.135.34; PORT=7777; DEST=data/BioMolDB_20261014; else REMOTE=psk6950@147.47.218.30; PORT=7994; DEST=/public_data02/BioMolDB_20261014; fi
  export REMOTE DEST SSH="ssh -o BatchMode=yes -p $PORT -i $HOME/.ssh/id_ed25519 -c aes128-gcm@openssh.com"
  ssh -n -o BatchMode=yes -p $PORT -i $HOME/.ssh/id_ed25519 $REMOTE "mkdir -p $DEST/lmdb/design/db_complex_pdb"
  tr '\t\n' '\0\0' < $HERE/units_${T}_rank0.txt | xargs -0 -n3 -P6 $HERE/send_one.sh
  $HERE/verify2.sh $T
done
