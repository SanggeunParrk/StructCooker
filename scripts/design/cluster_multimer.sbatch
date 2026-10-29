#!/bin/bash
#SBATCH -J design_fsm
#SBATCH -c 112
#SBATCH --mem=480g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 24:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/design/fsm_%j.out
# usage: sbatch cluster_multimer.sbatch DIR   (DIR/part_*.tar of multi-chain PDB files -> DIR/fs/clusters.tsv: representative, member, one line per complex)
cd "${SLURM_SUBMIT_DIR:?}"
D=$1; W=$D/fs; rm -rf $W; mkdir -p $W
s=$(date +%s)
# easy-multimercluster reads only its first input, so the tars are merged into one first
cp $(ls $D/part_*.tar | head -1) $W/all.tar
for t in $(ls $D/part_*.tar | tail -n +2); do tar --concatenate --file=$W/all.tar $t; done
echo "merged $(tar tf $W/all.tar | wc -l) structures in $(( $(date +%s) - s ))s"; s=$(date +%s)
pixi run foldseek easy-multimercluster $W/all.tar $W/res $W/tmp --threads 112 -c 0.9 -v 1
echo "easy-multimercluster $(( $(date +%s) - s ))s"
mv $W/res_cluster.tsv $W/clusters.tsv
rm -f $W/all.tar
echo "clusters $(cut -f1 $W/clusters.tsv | sort -u | wc -l) of $(wc -l < $W/clusters.tsv)"
