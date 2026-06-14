#!/bin/bash
#SBATCH -J easy_cluster_filter
#SBATCH --ntasks-per-node=1
#SBATCH -c 80
#SBATCH --mem=488g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/easy_cluster_filter.out
#SBATCH -e ./logs/easy_cluster_filter.err

set -euo pipefail

# =========================
# Usage
# =========================
# sbatch submits/easy_cluster_filter.sh <query.fasta> <key.fasta>
if [ $# -lt 2 ]; then
    echo "Usage: $0 <query.fasta> <key.fasta>" >&2
    exit 1
fi

QUERY=$(readlink -f "$1")
KEY=$(readlink -f "$2")
THREADS=80
MIN_LEN=200
MAX_LEN=640

# =========================
# working directory = query가 있는 디렉토리
# =========================
WD=$(dirname "$QUERY")
cd "$WD"

QUERY_NAME=$(basename "$QUERY")

# =========================
# 1) QUERY를 30% identity로 easy-cluster (FASTA -> FASTA)
#    출력:
#      clust_rep_seq.fasta   (대표 서열)
#      clust_all_seqs.fasta  (모든 서열, cluster 정보 포함)
#      clust_cluster.tsv     (rep \t member)
# =========================
mmseqs easy-cluster "$QUERY_NAME" clust cluster_tmp \
    --min-seq-id 0.3 \
    -c 0.8 \
    --cov-mode 0 \
    --threads "$THREADS"

# =========================
# 2) 대표 서열에 길이 필터 (MIN_LEN <= seq_len <= MAX_LEN)
# =========================
awk -v min="$MIN_LEN" -v max="$MAX_LEN" '
    function emit() {
        if (header != "" && length(seq) >= min && length(seq) <= max) {
            print header
            print seq
        }
    }
    /^>/ { emit(); header = $0; seq = ""; next }
    { seq = seq $0 }
    END   { emit() }
' clust_rep_seq.fasta > query_rep_len_filtered.fasta

# =========================
# 3) 검색: length-filtered 대표 서열 vs KEY
# =========================
mmseqs easy-search query_rep_len_filtered.fasta "$KEY" result.tsv search_tmp \
    -a --min-seq-id 0.3 \
    --threads "$THREADS"

# =========================
# 4) identity > 0.30인 query id 추출
# =========================
awk '$3 > 0.30 {print $1}' result.tsv | sort -u > remove_ids.txt

# =========================
# 5) 최종 fasta: 길이 필터 통과 + identity hit 제거
# =========================
awk 'NR==FNR {rm[$1]=1; next}
     /^>/ {
         id=$1
         sub(/^>/,"",id)
         keep = !(id in rm)
     }
     keep {print}' remove_ids.txt query_rep_len_filtered.fasta > query_final.fasta

# =========================
# 6) 요약
# =========================
TOTAL_QUERY=$(grep -c '^>' "$QUERY_NAME" || true)
REP_COUNT=$(grep -c '^>' clust_rep_seq.fasta || true)
REP_LEN_FILTERED=$(grep -c '^>' query_rep_len_filtered.fasta || true)
REMOVED=$(wc -l < remove_ids.txt || true)
FINAL=$(grep -c '^>' query_final.fasta || true)

echo "Done."
echo "Working directory : $WD"
echo "Query fasta       : $QUERY"
echo "Key fasta         : $KEY"
echo "Length filter     : ${MIN_LEN} <= len <= ${MAX_LEN}"
echo "Total query seqs  : $TOTAL_QUERY"
echo "Cluster reps      : $REP_COUNT"
echo "After length flt  : $REP_LEN_FILTERED"
echo "Removed by hit    : $REMOVED"
echo "Final seqs        : $FINAL"
echo "Output            : $WD/query_final.fasta"
