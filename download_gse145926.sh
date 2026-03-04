#!/usr/bin/env bash
# ============================================================
#  download_gse145926.sh
#  GSE145926 全サンプルのGEX(H5) + TCR(CSV) をダウンロード
#
#  Usage:
#    bash download_gse145926.sh            # 全サンプル
#    bash download_gse145926.sh C142       # 指定サンプルのみ
#
#  Run: conda activate alew の必要なし（bashのみ）
#
#  注意: ファイルは各GSMの個別サプリメントページから取得
#        /geo/samples/GSMxxxxnnn/GSMxxxx/suppl/
# ============================================================

set -uo pipefail

BASE_DIR="/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database"
DATA_DIR="${BASE_DIR}/data"
FTP_SAMPLES_BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples"

# ── サンプル定義 ─────────────────────────────────────────────
# 形式: "SAMPLE_ID|GSM_GEX|GSM_TCR"
# TCRなしのサンプルはGSM_TCRを"-"にする
declare -a SAMPLES=(
    "C141|GSM4339769|GSM4385990"
    "C142|GSM4339770|GSM4385991"
    "C143|GSM4339771|GSM4385992"
    "C144|GSM4339772|GSM4385993"
    "C145|GSM4339773|GSM4385994"
    "C146|GSM4339774|GSM4385995"
    "C51|GSM4475048|-"
    "C52|GSM4475049|-"
    "C100|GSM4475050|-"
    "C148|GSM4475051|GSM4475054"
    "C149|GSM4475052|GSM4475055"
    "C152|GSM4475053|GSM4475056"
)

# ── 対象サンプルのフィルタ ────────────────────────────────────
TARGET="${1:-}"   # 引数があれば指定サンプルのみ

# ── GSMのサブディレクトリ名を計算 ────────────────────────────
# 例: GSM4339769 → GSM4339nnn
gsm_to_dir() {
    local gsm="$1"
    # 下3桁を nnn に置換
    echo "${gsm:0:-3}nnn"
}

# ── ダウンロード関数 ──────────────────────────────────────────
download_file() {
    local url="$1"
    local dest="$2"
    if [[ -f "$dest" && -s "$dest" ]]; then
        echo "  [SKIP] Already exists: $(basename "$dest")"
        return 0
    fi
    echo "  [DL]   $(basename "$dest") from $(dirname "$url")"
    rm -f "${dest}.tmp"
    if wget -q --show-progress --timeout=300 --tries=3 \
            -O "${dest}.tmp" "$url"; then
        mv "${dest}.tmp" "$dest"
        echo "  [OK]   $(basename "$dest") ($(du -sh "$dest" | cut -f1))"
    else
        rm -f "${dest}.tmp"
        echo "  [FAIL] Failed to download: $(basename "$dest")"
        return 1
    fi
}

# ── メイン処理 ────────────────────────────────────────────────
echo "================================================"
echo " GSE145926 Sample Downloader"
echo " Target: ${TARGET:-ALL}"
echo "================================================"

PROCESSED=0
FAILED=0

for entry in "${SAMPLES[@]}"; do
    IFS='|' read -r sample_id gsm_gex gsm_tcr <<< "$entry"
    full_id="GSE145926_${sample_id}"

    # 引数で絞り込み
    if [[ -n "$TARGET" && "$sample_id" != "$TARGET" ]]; then
        continue
    fi

    echo ""
    echo ">>> Sample: ${sample_id} (${full_id})"
    out_dir="${DATA_DIR}/${full_id}"
    mkdir -p "$out_dir"

    SAMPLE_FAILED=0

    # ── GEX H5 ──────────────────────────────────────────────
    gsm_gex_dir=$(gsm_to_dir "$gsm_gex")
    h5_filename="${gsm_gex}_${sample_id}_filtered_feature_bc_matrix.h5"
    h5_url="${FTP_SAMPLES_BASE}/${gsm_gex_dir}/${gsm_gex}/suppl/${h5_filename}"
    h5_dest="${out_dir}/filtered_feature_bc_matrix.h5"
    if ! download_file "$h5_url" "$h5_dest"; then
        ((FAILED++)) || true
        SAMPLE_FAILED=1
    fi

    # ── TCR CSV ─────────────────────────────────────────────
    if [[ "$gsm_tcr" != "-" ]]; then
        gsm_tcr_dir=$(gsm_to_dir "$gsm_tcr")
        tcr_gz_filename="${gsm_tcr}_${sample_id}_filtered_contig_annotations.csv.gz"
        tcr_gz_url="${FTP_SAMPLES_BASE}/${gsm_tcr_dir}/${gsm_tcr}/suppl/${tcr_gz_filename}"
        tcr_gz_dest="${out_dir}/filtered_contig_annotations_tcr.csv.gz"
        tcr_dest="${out_dir}/filtered_contig_annotations_tcr.csv"
        if download_file "$tcr_gz_url" "$tcr_gz_dest"; then
            if [[ -f "$tcr_gz_dest" && ! -f "$tcr_dest" ]]; then
                echo "  [GZ]   Decompressing TCR CSV..."
                gunzip -kf "$tcr_gz_dest"
                # gunzip は .gz を除いたファイルを同ディレクトリに作成
                mv "${out_dir}/${tcr_gz_filename%.gz}" "$tcr_dest" 2>/dev/null || true
            fi
        else
            ((FAILED++)) || true
            SAMPLE_FAILED=1
        fi
    else
        echo "  [INFO] No TCR data for ${sample_id}"
    fi

    if [[ "$SAMPLE_FAILED" == "0" ]]; then
        echo "  [OK]   ${full_id} ready in: ${out_dir}"
        ((PROCESSED++)) || true
    fi
done

echo ""
echo "================================================"
echo " Download complete: ${PROCESSED} sample(s) processed"
if [[ "$FAILED" -gt 0 ]]; then
    echo " Failed: ${FAILED} file(s)"
fi
echo " Data directory   : ${DATA_DIR}"
echo "================================================"
echo ""
echo "Next step: Run batch processing"
echo "  bash batch_process_gse145926.sh"
