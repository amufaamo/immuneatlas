#!/usr/bin/env bash
# ============================================================
#  batch_process_gse145926.sh
#  GSE145926 全サンプルを順次 Rscript で処理する
#
#  Usage:
#    bash batch_process_gse145926.sh          # 全サンプル
#    bash batch_process_gse145926.sh C142     # 指定サンプルのみ
#
#  前提: conda activate alew は自動で行う (conda run 経由)
# ============================================================

set -euo pipefail

BASE_DIR="/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database"
RSCRIPT="${BASE_DIR}/process_geo_sample.R"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

# 処理対象サンプル (TCRありなし関係なく全て)
ALL_SAMPLES=(C141 C142 C143 C144 C145 C146 C51 C52 C100 C148 C149 C152)

TARGET="${1:-}"  # 引数があれば指定サンプルのみ

echo "======================================================="
echo " GSE145926 Batch Processor"
echo " Target  : ${TARGET:-ALL (${#ALL_SAMPLES[@]} samples)}"
echo " R script: ${RSCRIPT}"
echo " Logs    : ${LOG_DIR}"
echo "======================================================="

DONE=0
FAILED=0
SKIPPED_LIST=()

for sample_id in "${ALL_SAMPLES[@]}"; do
    # 引数フィルタ
    if [[ -n "$TARGET" && "$sample_id" != "$TARGET" ]]; then
        continue
    fi

    full_id="GSE145926_${sample_id}"
    data_dir="${BASE_DIR}/data/${full_id}"
    out_json="${BASE_DIR}/output/umap_data_${full_id}.json"
    log_file="${LOG_DIR}/process_${sample_id}.log"

    echo ""
    echo ">>> [$(date '+%H:%M:%S')] Processing ${sample_id}..."

    # データが存在しない場合はスキップ
    if [[ ! -d "$data_dir" || ! -f "${data_dir}/filtered_feature_bc_matrix.h5" ]]; then
        echo "  [SKIP] Data not found: ${data_dir}"
        echo "         Run: bash download_gse145926.sh ${sample_id}"
        SKIPPED_LIST+=("$sample_id")
        continue
    fi

    # 既に処理済みの場合はスキップ (--force フラグがなければ)
    if [[ -f "$out_json" && "${FORCE:-0}" != "1" ]]; then
        echo "  [SKIP] Already processed (output exists): ${out_json}"
        echo "         To re-run: FORCE=1 bash batch_process_gse145926.sh ${sample_id}"
        SKIPPED_LIST+=("${sample_id}(cached)")
        continue
    fi

    # Rscript 実行 (conda run 経由)
    echo "  [RUN]  conda run -n alew Rscript ... --sample ${sample_id}"
    if conda run -n alew Rscript "$RSCRIPT" --sample "$sample_id" \
            > "$log_file" 2>&1; then
        echo "  [OK]   Done. Log: ${log_file}"
        ((DONE++)) || true
    else
        echo "  [FAIL] Error in ${sample_id}! Check log: ${log_file}"
        tail -20 "$log_file" | sed 's/^/         /'
        ((FAILED++)) || true
    fi
done

echo ""
echo "======================================================="
echo " Batch complete"
echo "   Processed : ${DONE}"
echo "   Failed    : ${FAILED}"
if [[ ${#SKIPPED_LIST[@]} -gt 0 ]]; then
    echo "   Skipped   : ${SKIPPED_LIST[*]}"
fi
echo "======================================================="

if [[ "$DONE" -gt 0 ]]; then
    echo ""
    echo "Next step: Regenerate studies.json"
    echo "  bash generate_studies_json.sh"
fi
