#!/bin/bash
#SBATCH -J design_fs
#SBATCH -c 112
#SBATCH --mem=480g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 12:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/design/fs_%j.out
# usage: sbatch cluster_tar.sbatch DIR   (DIR/part_*.tar from build_stats -> DIR/fs/clusters.tsv, foldseek defaults, -c 0.9)
cd "${SLURM_SUBMIT_DIR:?}"
D=$1; W=$D/fs; mkdir -p $W; FS="pixi run foldseek"
s=$(date +%s)
$FS createdb $D/part_*.tar $W/db --threads 112 -v 1
echo "createdb $(( $(date +%s) - s ))s"; s=$(date +%s)
rm -rf $W/tmp; $FS cluster $W/db $W/clu $W/tmp --threads 112 -c 0.9 --split-memory-limit 400G -v 1
echo "cluster $(( $(date +%s) - s ))s"
$FS createtsv $W/db $W/db $W/clu $W/clusters.tsv -v 1
echo "clusters $(cut -f1 $W/clusters.tsv | sort -u | wc -l) of $(wc -l < $W/clusters.tsv)"
