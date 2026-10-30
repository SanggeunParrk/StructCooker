#!/bin/bash
#SBATCH -J tfix_p12
#SBATCH -c 32
#SBATCH --mem=64g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 2:00:00
#SBATCH --exclude=node01,node04,node05,node07
#SBATCH -o logs/template_fix/p12_%j.out
cd "${SLURM_SUBMIT_DIR:?}"
export BIOMOL_ROOT=/data/psk6950/template_fix/root
S=.pixi/envs/default/bin/structcooker
for s in template/template_candidates template/template_seqs metadata/seqid_template_filelist; do echo "== $s"; $S build $s || exit 1; done
wc -l $BIOMOL_ROOT/materials/intermediate/worklists/*.tsv $BIOMOL_ROOT/materials/intermediate/worklists/*.txt
