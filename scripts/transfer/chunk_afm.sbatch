#!/bin/bash
#SBATCH -J afm_chunks
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -c 8
#SBATCH --mem=16g
#SBATCH --time=1-00:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o /home/psk6950/StructCooker/logs/transfer/afm_chunks_%A_%a.out
HERE=/home/psk6950/StructCooker/logs/transfer
echo "rank $SLURM_ARRAY_TASK_ID $(hostname) $(wc -l < $HERE/chunks_afm_rank$SLURM_ARRAY_TASK_ID.txt) chunks $(date +%T)"
tr '\t\n' '\0\0' < $HERE/chunks_afm_rank$SLURM_ARRAY_TASK_ID.txt | xargs -0 -n3 -P14 $HERE/chunk_one.sh
echo "done $(date +%T)"
