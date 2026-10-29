#!/bin/bash
#SBATCH -J if_cluster_sim
#SBATCH -c 56
#SBATCH --mem=200g
#SBATCH -p cpu-long
#SBATCH --qos=cpu-long-q
#SBATCH -t 12:00:00
#SBATCH --exclude=node01,node04,node07
#SBATCH -o logs/analysis/if_cluster_sim_%j.out
# usage: sbatch interface_cluster_similarity.sbatch OUT [--n-clusters N] [--max-members M]
cd "${SLURM_SUBMIT_DIR:?}"
OUT=$1; shift
pixi run python -I scripts/analysis/interface_cluster_similarity.py \
  --interfaces /data/psk6950/train_item/pdb \
  --db /data/shared/cssb_data/BioMol_clean/lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb \
  --out "$OUT" --workers 56 "$@" 2>&1 | grep -v WARN
