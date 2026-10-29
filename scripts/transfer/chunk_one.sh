#!/bin/bash
# chunk_one.sh REL OFFSET LEN -- write bytes [OFFSET, OFFSET+LEN) of the source file into the same file on H100 (in place, no truncation)
REL=$1; OFF=$2; LEN=$3
SRC=/data/shared/cssb_data/BioMol_clean/lmdb/afdb_multimer/$REL
DST=data/BioMolDB_20261014/lmdb/afdb_multimer/$REL
DONE=/data/psk6950/transfer/chunk_done; mkdir -p $DONE
M=$DONE/$(echo "$REL" | tr '/' '_').$OFF
[ -e "$M" ] && exit 0
sleep $((RANDOM % 8))
for t in 1 2 3 4; do
  if dd if=$SRC bs=16M iflag=skip_bytes,count_bytes skip=$OFF count=$LEN status=none | \
     ssh -o BatchMode=yes -p 7777 -i $HOME/.ssh/id_ed25519 -c aes128-gcm@openssh.com psk6950@147.47.135.34 "dd of=$DST bs=16M oflag=seek_bytes seek=$OFF conv=notrunc status=none"; then
    [ "${PIPESTATUS[0]}" = 0 ] && touch "$M" && echo "ok $REL $OFF $(date +%T)" && exit 0
  fi
  sleep 15
done
echo "FAILED $REL $OFF"; exit 1
