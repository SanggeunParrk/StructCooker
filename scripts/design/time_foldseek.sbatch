#!/bin/bash
#SBATCH -J fs_time
#SBATCH -c 32
#SBATCH --mem=64g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 2:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/design/fs_time_%j.out
# usage: sbatch time_foldseek.sbatch PDB_DIR WORK_DIR   (times createdb and several cluster settings)
cd "${SLURM_SUBMIT_DIR:?}"
P=$1; W=$2; mkdir -p $W; FS="pixi run foldseek"
t() { local s=$(date +%s); "$@" >$W/last.log 2>&1; echo "  rc=$? $(( $(date +%s) - s ))s"; }
echo createdb; t $FS createdb $P $W/db --threads 32
echo "cluster default (3Di prefilter, local alignment 3Di+AA, c0.9)"; rm -rf $W/tmp1; t $FS cluster $W/db $W/c1 $W/tmp1 --threads 32 -c 0.9 -v 1
$FS createtsv $W/db $W/db $W/c1 $W/c1.tsv >/dev/null 2>&1; echo "  clusters $(cut -f1 $W/c1.tsv | sort -u | wc -l)"
echo "cluster tmalign (alignment-type 1)"; rm -rf $W/tmp2; t $FS cluster $W/db $W/c2 $W/tmp2 --threads 32 -c 0.9 --alignment-type 1 --tmscore-threshold 0.5 -v 1
$FS createtsv $W/db $W/db $W/c2 $W/c2.tsv >/dev/null 2>&1; echo "  clusters $(cut -f1 $W/c2.tsv | sort -u | wc -l)"
echo "cluster default, c0.8 cov-mode 0 with tm filter"; rm -rf $W/tmp3; t $FS cluster $W/db $W/c3 $W/tmp3 --threads 32 -c 0.8 --alignment-type 2 --tmscore-threshold 0.5 -v 1
$FS createtsv $W/db $W/db $W/c3 $W/c3.tsv >/dev/null 2>&1; echo "  clusters $(cut -f1 $W/c3.tsv | sort -u | wc -l)"
