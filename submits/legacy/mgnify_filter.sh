#!/bin/bash
#SBATCH -J mgnify_filter
#SBATCH --ntasks-per-node=1
#SBATCH -c 80
#SBATCH --mem=488g
#SBATCH -p cpu
#SBATCH -w node02
#SBATCH -o ./logs/mgnify_filter.out
#SBATCH -e ./logs/mgnify_filter.err

set -euo pipefail

# =========================
# 설정
# =========================
WD=/public_data02/distillation_structures/mgnify
QUERY=/public_data02/distillation_structures/mgnify/mgy_clusters_2022_05.fa
KEY=/public_data02/BioMolDB_20260224/fasta/merged.fasta
THREADS=80
MIN_LEN=200
MAX_LEN=640

# =========================
# working directory 준비
# =========================
mkdir -p "$WD"
cd "$WD"

QUERY_NAME=$(basename "$QUERY")
KEY_NAME=$(basename "$KEY")

# =========================
# 입력 파일 준비
# - QUERY는 WD 밖에 있으므로 복사
# - KEY는 이미 WD 안에 있을 수 있으므로 필요할 때만 복사
# =========================
if [ "$(readlink -f "$QUERY")" != "$(readlink -f "$WD/$QUERY_NAME")" ]; then
    cp "$QUERY" "$WD/$QUERY_NAME"
fi

if [ "$(readlink -f "$KEY")" != "$(readlink -f "$WD/$KEY_NAME")" ]; then
    cp "$KEY" "$WD/$KEY_NAME"
fi

# =========================
# 1) FASTA -> MMseqs DB
# =========================
mmseqs createdb "$QUERY_NAME" queryDB
mmseqs createdb "$KEY_NAME"   keyDB

# =========================
# 1-1) QUERY를 30% identity로 클러스터링 → 대표서열 DB
# clusterDB.index의 key(첫 컬럼)가 각 cluster의 대표 key
# =========================
mmseqs cluster queryDB clusterDB cluster_tmp \
    --min-seq-id 0.3 \
    -c 0.8 \
    --cov-mode 0 \
    -s 6 \
    --threads "$THREADS"

awk '{print $1}' clusterDB.index > rep_keys.txt
mmseqs createsubdb rep_keys.txt queryDB queryDB_rep

# =========================
# 1-2) 대표서열 DB에 길이 필터 (MIN_LEN <= seq_len <= MAX_LEN)
# mmseqs .index의 length는 seq 뒤에 '\n\0'이 붙어 있어 실제 길이보다 +2
# =========================
awk -v min="$MIN_LEN" -v max="$MAX_LEN" \
    '($3 - 2) >= min && ($3 - 2) <= max {print $1}' \
    queryDB_rep.index > keep_keys.txt

mmseqs createsubdb keep_keys.txt queryDB_rep queryDB_rep_lf

# =========================
# 2) 검색 (clustered + length-filtered query vs key)
# -a : alignment 계산
# --min-seq-id 0.3 : 30% 이상 hit 유지
# =========================
mmseqs search queryDB_rep_lf keyDB resultDB tmp \
    -a --min-seq-id 0.3 -s 6 \
    --threads "$THREADS"

# =========================
# 3) 결과를 TSV로 변환
# =========================
mmseqs convertalis queryDB_rep_lf keyDB resultDB result.tsv

# =========================
# 4) 제거할 query ID 추출 (identity > 0.30)
# =========================
awk '$3 > 0.30 {print $1}' result.tsv | sort -u > remove_ids.txt

# =========================
# 5) 클러스터 대표 + 길이 필터 통과한 fasta 추출 후 identity hit 제거
# FASTA header 첫 토큰 기준
# =========================
mmseqs convert2fasta queryDB_rep_lf query_rep_len_filtered.fasta

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
REP_COUNT=$(wc -l < rep_keys.txt || true)
REP_LEN_FILTERED=$(grep -c '^>' query_rep_len_filtered.fasta || true)
REMOVED=$(wc -l < remove_ids.txt || true)
FINAL=$(grep -c '^>' query_final.fasta || true)

echo "Done."
echo "Working directory : $WD"
echo "Query fasta       : $QUERY_NAME"
echo "Key fasta         : $KEY_NAME"
echo "Length filter     : ${MIN_LEN} <= len <= ${MAX_LEN}"
echo "Total query seqs  : $TOTAL_QUERY"
echo "Cluster reps      : $REP_COUNT"
echo "After length flt  : $REP_LEN_FILTERED"
echo "Removed by hit    : $REMOVED"
echo "Final seqs        : $FINAL"
echo "Output            : $WD/query_final.fasta"
