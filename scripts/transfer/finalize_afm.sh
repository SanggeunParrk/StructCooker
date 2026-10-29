#!/bin/bash
# size check + 24 random 64 MiB regions per file (md5 local vs H100) for the three chunk-written AFDB files
SRC=/data/shared/cssb_data/BioMol_clean/lmdb/afdb_multimer
R="ssh -o BatchMode=yes -p 7777 -i $HOME/.ssh/id_ed25519 psk6950@147.47.135.34"
for rel in msa/msa_afm_d2k.lmdb/data.mdb template/template_afm.lmdb/data.mdb msa/msa_afm.lmdb/data.mdb; do
  ls=$(stat -c %s $SRC/$rel); rs=$($R "stat -c %s data/BioMolDB_20261014/lmdb/afdb_multimer/$rel")
  [ "$ls" = "$rs" ] && echo "size ok $rel $ls" || { echo "SIZE MISMATCH $rel $ls $rs"; continue; }
  bad=0
  for i in $(seq 1 24); do
    off=$(( (RANDOM * 32768 + RANDOM) % (ls / 67108864) ))
    l=$(dd if=$SRC/$rel bs=64M skip=$off count=1 status=none | md5sum | cut -c1-12)
    r=$($R "dd if=data/BioMolDB_20261014/lmdb/afdb_multimer/$rel bs=64M skip=$off count=1 status=none | md5sum" | cut -c1-12)
    [ "$l" = "$r" ] || { echo "REGION MISMATCH $rel block $off"; bad=1; }
  done
  [ $bad = 0 ] && echo "regions ok $rel"
  $R "touch -d @$(stat -c %Y $SRC/$rel) data/BioMolDB_20261014/lmdb/afdb_multimer/$rel"
done
