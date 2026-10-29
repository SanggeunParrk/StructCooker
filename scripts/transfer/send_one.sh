#!/bin/bash
# send_one.sh ROOT REL PREFIX   (env: SSH REMOTE DEST) -- one file, relative path kept under DEST/PREFIX, 3 tries
ROOT=$1; REL=$2; PRE=$3
sleep $((RANDOM % 6)); cd "$ROOT" || exit 1
for t in 1 2 3; do
  if rsync -aWR --partial --exclude=lock.mdb -e "$SSH" "./$REL" "$REMOTE:$DEST/$PRE/"; then echo "sent $PRE/$REL $(date +%T)"; exit 0; fi
  sleep 20
done
echo "FAILED $PRE/$REL"; exit 1
