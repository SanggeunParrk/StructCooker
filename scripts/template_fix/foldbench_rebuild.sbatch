#!/bin/bash
#SBATCH -J fb_rebuild
#SBATCH -c 8
#SBATCH --mem=64g
#SBATCH -p cpu-short
#SBATCH --qos=cpu-short-q
#SBATCH -t 1:00:00
#SBATCH --exclude=node01,node04,node05,node07
#SBATCH -o logs/template_fix/fb_rebuild_%j.out
cd "${SLURM_SUBMIT_DIR:?}"
R=/data/shared/cssb_data/BioMol_clean
O=/data/psk6950/foldbench_rebuild_20261029
rm -rf $O
S=.pixi/envs/default/bin/structcooker
$S foldbench build $R/materials/reference/foldbench $O --release $R --ccd $R/lmdb/ccd/ccd_inference.lmdb --fallback $R/lmdb/foldbench/msa_edge.lmdb --local && \
$S foldbench verify $R/materials/reference/foldbench $O
