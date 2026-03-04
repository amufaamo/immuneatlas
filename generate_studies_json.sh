#!/usr/bin/env bash
# ============================================================
#  generate_studies_json.sh
#  output/ 以下の umap_data_*.json から studies.json を自動生成
#  Run after batch_process_gse145926.sh
# ============================================================

BASE_DIR="/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database"
OUTPUT_DIR="${BASE_DIR}/output"
WEBAPP_DIR="${BASE_DIR}/webapp"
OUT_FILE="${OUTPUT_DIR}/studies.json"

# ── サンプルメタデータ定義 ───────────────────────────────────
# "sample_key|display_title|disease|severity|assay|organ|data_url_suffix|viewer_param"
declare -A META_TITLE META_DISEASE META_SEVERITY META_ASSAY META_ORGAN META_GROUP META_GROUP_TITLE

# ALEW local sample
META_TITLE["mvp_sample"]="ALEW Cohort — MVP Single-cell Immune Repertoire"
META_DISEASE["mvp_sample"]="Diseased"
META_SEVERITY["mvp_sample"]=""
META_ASSAY["mvp_sample"]="10x 5' GEX + TCR/BCR"
META_ORGAN["mvp_sample"]="PBMC"
META_GROUP["mvp_sample"]="mvp_sample"
META_GROUP_TITLE["mvp_sample"]="ALEW Cohort — MVP Single-cell Immune Repertoire"

# GSE145926 samples (Liao et al. Nature Medicine 2020)
GSE_TITLE="Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19"
GSE_DESC="Respiratory immune characteristics associated with COVID-19 severity. Bronchoalveolar lavage fluid (BALF) profiled by scRNA-seq. Liao M et al. Nature Medicine 2020. DOI:10.1038/s41591-020-0901-9 (GEO: GSE145926)"

for sid in C141 C142 C143 C144 C145 C146 C51 C52 C100 C148 C149 C152; do
    key="gse145926_${sid,,}"
    META_TITLE[$key]="${GSE_TITLE} [${sid}]"
    META_DISEASE[$key]="COVID-19"
    META_ASSAY[$key]="10x 5' GEX + TCR"
    META_ORGAN[$key]="BALF"
    META_GROUP[$key]="gse145926"
    META_GROUP_TITLE[$key]="${GSE_TITLE} (GSE145926)"
done
# Severity annotation (Liao et al.)
META_SEVERITY["gse145926_c141"]="Moderate"
META_SEVERITY["gse145926_c142"]="Moderate"
META_SEVERITY["gse145926_c143"]="Moderate"
META_SEVERITY["gse145926_c144"]="Severe"
META_SEVERITY["gse145926_c145"]="Severe"
META_SEVERITY["gse145926_c146"]="Severe"
META_SEVERITY["gse145926_c51"]="Severe"
META_SEVERITY["gse145926_c52"]="Moderate"
META_SEVERITY["gse145926_c100"]="Moderate"
META_SEVERITY["gse145926_c148"]="Severe"
META_SEVERITY["gse145926_c149"]="Severe"
META_SEVERITY["gse145926_c152"]="Severe"

# ── JSON 生成 ────────────────────────────────────────────────
echo "Scanning output files in: ${OUTPUT_DIR}"

# JSON 開始
cat > "$OUT_FILE" << 'JSONSTART'
[
JSONSTART

FIRST=true

write_entry() {
    local key="$1"
    local data_file="$2"
    local viewer_param="$3"
    local is_merged="${4:-false}"
    local qs2_file=""
    
    if [[ "$is_merged" == "true" ]]; then
        qs2_file="${data_file%.json}.qs2"
        # qs2_file が umap_data_gse145926_merged.qs2 みたいな形になるが、
        # 実際は seurat_gse145926_merged.qs2 かもしれないのでチェック
        local dir=$(dirname "$data_file")
        local base=$(basename "$data_file")
        local s_qs2="${dir}/seurat_${base#umap_data_}"
        s_qs2="${s_qs2%.json}.qs2"
        if [[ -f "$s_qs2" ]]; then qs2_file="$s_qs2"; fi
    else
        qs2_file="${OUTPUT_DIR}/seurat_$(basename "${data_file%.json}" | sed 's/umap_data_//').qs2"
    fi

    local title="${META_TITLE[$key]:-Unknown Sample ($key)}"
    local disease="${META_DISEASE[$key]:-Unknown}"
    local severity="${META_SEVERITY[$key]:-}"
    local assay="${META_ASSAY[$key]:-scRNA-seq}"
    local organ="${META_ORGAN[$key]:-Unknown}"
    local group="${META_GROUP[$key]:-$key}"
    local group_title="${META_GROUP_TITLE[$key]:-$title}"
    local qs2_path=""
    if [[ -f "$qs2_file" ]]; then
        qs2_path="/output/$(basename "$qs2_file")"
    fi

    # Diversity データの取得 (jq or python)
    local div_info=""
    local div_file="${OUTPUT_DIR}/diversity_summary.json"
    if [[ -f "$div_file" ]]; then
        # key にマッチするエントリを取得 (GSE145926_C141 形式か mvp_sample 以前のキーかチェック)
        local search_key="$key"
        if [[ "$key" == "mvp_sample" ]]; then search_key="MVP_Sample"; fi
        if [[ "$key" == "gse145926_merged" ]]; then search_key="gse145926_merged"; fi
        if [[ "$key" == gse145926_c* ]]; then
            local sid="${key#gse145926_}"
            sid="${sid^^}"
            search_key="GSE145926_${sid}"
        fi
        
        div_info=$(python3 -c "import sys,json; d=json.load(open('$div_file')); k='$search_key'; print(json.dumps(d.get(k, {})))")
    fi

    if [[ "$FIRST" != "true" ]]; then
        echo "," >> "$OUT_FILE"
    fi
    FIRST=false

    echo "  {" >> "$OUT_FILE"
    echo "    \"id\": \"${key}\"," >> "$OUT_FILE"
    echo "    \"title\": $(echo "$title" | python3 -c "import sys,json; print(json.dumps(sys.stdin.read().strip()))")," >> "$OUT_FILE"
    echo "    \"description\": $(echo "${GSE_DESC:-}" | python3 -c "import sys,json; print(json.dumps(sys.stdin.read().strip()))")," >> "$OUT_FILE"
    echo "    \"disease\": \"${disease}\"," >> "$OUT_FILE"
    echo "    \"severity\": \"${severity}\"," >> "$OUT_FILE"
    echo "    \"species\": \"Homo sapiens\"," >> "$OUT_FILE"
    echo "    \"assay\": \"${assay}\"," >> "$OUT_FILE"
    echo "    \"organ\": \"${organ}\"," >> "$OUT_FILE"
    echo "    \"studyGroup\": \"${group}\"," >> "$OUT_FILE"
    echo "    \"studyGroupTitle\": $(echo "$group_title" | python3 -c "import sys,json; print(json.dumps(sys.stdin.read().strip()))")," >> "$OUT_FILE"
    
    if [[ "$is_merged" == "true" ]]; then
        echo "    \"isMerged\": true," >> "$OUT_FILE"
        echo "    \"mergedStudyId\": \"${viewer_param}\"," >> "$OUT_FILE"
        echo "    \"viewerUrl\": \"/webapp/viewer.html?study=${viewer_param}\"," >> "$OUT_FILE"
    else
        echo "    \"viewerUrl\": \"/webapp/viewer.html?sample=${viewer_param}\"," >> "$OUT_FILE"
    fi
    
    # Diversity info integration
    if [[ -n "$div_info" && "$div_info" != "{}" ]]; then
        local n_cells=$(echo "$div_info" | python3 -c "import sys,json; v=json.load(sys.stdin).get('n_cells'); print(json.dumps(v))")
        local n_tcr=$(echo "$div_info" | python3 -c "import sys,json; v=json.load(sys.stdin).get('n_tcr'); print(json.dumps(v))")
        local shannon=$(echo "$div_info" | python3 -c "import sys,json; v=json.load(sys.stdin).get('shannon'); print(json.dumps(v))")
        
        echo "    \"cellCount\": ${n_cells:-0}," >> "$OUT_FILE"
        echo "    \"nClonotypes\": ${n_tcr:-0}," >> "$OUT_FILE"
        echo "    \"shannonIndex\": ${shannon:-null}," >> "$OUT_FILE"
    fi

    echo "    \"dataUrl\": \"/output/$(basename "$data_file")\"," >> "$OUT_FILE"
    echo "    \"qs2Url\": \"${qs2_path}\"" >> "$OUT_FILE"
    echo "  }" >> "$OUT_FILE"
}

# 1. MVP sample (固定)
MVP_JSON="${OUTPUT_DIR}/umap_data.json"
if [[ -f "$MVP_JSON" ]]; then
    write_entry "mvp_sample" "$MVP_JSON" "mvp"
fi

# 2. Merged Studies
MERGED_JSON="${OUTPUT_DIR}/umap_data_gse145926_merged.json"
if [[ -f "$MERGED_JSON" ]]; then
    # Merged 用のメタデータをセット (既存のグループに統合)
    META_TITLE["gse145926_merged"]="${GSE_TITLE} — Merged Study View"
    META_GROUP["gse145926_merged"]="gse145926"
    META_GROUP_TITLE["gse145926_merged"]="${GSE_TITLE} (GSE145926)"
    META_DISEASE["gse145926_merged"]="COVID-19 (Integrated)"
    META_ORGAN["gse145926_merged"]="BALF"
    META_ASSAY["gse145926_merged"]="scRNA-seq (Harmony)"
    write_entry "gse145926_merged" "$MERGED_JSON" "gse145926" "true"
fi

# 3. GSE145926 individual samples
for sid in C141 C142 C143 C144 C145 C146 C51 C52 C100 C148 C149 C152; do
    full_id="GSE145926_${sid}"
    json_file="${OUTPUT_DIR}/umap_data_${full_id}.json"
    key="gse145926_${sid,,}"
    if [[ -f "$json_file" ]]; then
        write_entry "$key" "$json_file" "gse145926_${sid,,}"
        echo "  Added: ${full_id}"
    fi
done

echo "" >> "$OUT_FILE"
echo "]" >> "$OUT_FILE"

echo ""
echo "Generated: ${OUT_FILE}"
echo "Studies included: $(grep -c '"id"' "$OUT_FILE") entries"
