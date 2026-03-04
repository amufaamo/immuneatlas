#!/usr/bin/env bash
# ダウンロード→バッチ処理→studies.json生成 の一括パイプライン
set -uo pipefail

BASE_DIR="/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "======================================"
echo "[$(date)] Pipeline START"
echo "======================================"

echo "[$(date)] Step 1: Downloading all samples..."
bash "${BASE_DIR}/download_gse145926.sh" 2>&1 | tee "${LOG_DIR}/dl_all.log"

echo ""
echo "[$(date)] Step 2: Batch processing with Rscript..."
bash "${BASE_DIR}/batch_process_gse145926.sh" 2>&1 | tee "${LOG_DIR}/batch_all.log"

echo ""
echo "[$(date)] Step 3: Generating studies.json..."
bash "${BASE_DIR}/generate_studies_json.sh" 2>&1 | tee "${LOG_DIR}/studies_json.log"

echo ""
echo "======================================"
echo "[$(date)] Pipeline COMPLETE"
echo "======================================"
