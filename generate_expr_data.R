# ============================================================
#  generate_expr_data.R
#  既存の seurat_*.qs2 から expr_data_*.json を生成する
#  Usage:
#    Rscript generate_expr_data.R --sample C141
#    Rscript generate_expr_data.R --merged gse145926
#    Rscript generate_expr_data.R --mvp
#    Rscript generate_expr_data.R --all   # 全個体サンプル（merged以外）一括処理
#  Run in: conda activate alew
# ============================================================

library(Seurat)
library(jsonlite)
library(qs2)

# ── コマンドライン引数パース ────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

sample_id <- NULL
merged_id <- NULL
process_all <- FALSE
process_mvp <- FALSE

for (i in seq_along(args)) {
    if (args[i] == "--sample" && i < length(args)) {
        sample_id <- args[i + 1]
    }
    if (args[i] == "--merged" && i < length(args)) {
        merged_id <- args[i + 1]
    }
    if (args[i] == "--all") {
        process_all <- TRUE
    }
    if (args[i] == "--mvp") {
        process_mvp <- TRUE
    }
}

BASE_DIR <- "/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database"
OUTPUT_DIR <- file.path(BASE_DIR, "output")
N_GENES <- 3000

# ── 1サンプル処理関数 ────────────────────────────────────
process_file <- function(qs2_path, out_json_path, label) {
    if (!file.exists(qs2_path)) {
        message("  [SKIP] qs2 not found: ", qs2_path)
        return(FALSE)
    }
    if (file.exists(out_json_path)) {
        message("  [SKIP] expr_data already exists: ", out_json_path)
        return(TRUE)
    }

    message("======================================================")
    message(" Processing: ", label)
    message(" Loading qs2: ", qs2_path)
    message("======================================================")

    seu <- qs_read(qs2_path)
    message("  Cells: ", ncol(seu))

    # Seurat v5: 複数レイヤーがある場合は JoinLayers する
    # (GetAssayData が 1 レイヤーしか受け付けないため)
    if (packageVersion("Seurat") >= "5.0.0") {
        # Assay を取得
        default_assay <- DefaultAssay(seu)
        # レイヤー名を確認
        lyrs <- Layers(seu[[default_assay]])
        if (length(grep("^data", lyrs)) > 1) {
            message("  Detected multiple data layers. Joining layers...")
            seu <- JoinLayers(seu)
        }
    }

    # Variable Features
    var_features <- VariableFeatures(seu)
    if (length(var_features) == 0) {
        message("  Running FindVariableFeatures...")
        seu <- FindVariableFeatures(seu, nfeatures = N_GENES)
        var_features <- VariableFeatures(seu)
    }
    top_genes <- head(var_features, N_GENES)
    message("  Genes to export: ", length(top_genes))

    # log-normalized マトリクス (data slot)
    expr_mat <- GetAssayData(seu, layer = "data")[top_genes, , drop = FALSE]
    cell_barcodes <- colnames(seu)

    # スパース辞書形式に変換 (高速化版)
    expr_list <- list()
    pb_total <- length(top_genes)

    for (i in seq_along(top_genes)) {
        gene <- top_genes[i]
        # expr_mat[i, ] はスパースベクトル
        row_data <- expr_mat[i, ]

        # 非ゼロのインデックス
        nonzero_idx <- which(row_data > 0)
        if (length(nonzero_idx) == 0) next

        # 値を取得して丸める
        vals_nonzero <- as.numeric(row_data[nonzero_idx])
        vals_nonzero <- round(vals_nonzero, 3)

        # スパース行リスト: [ [0-based indices], [values] ]
        # jsonliteが [ [1,2], [1.2,3.4] ] のように保存する
        expr_list[[gene]] <- list(as.integer(nonzero_idx - 1), vals_nonzero)

        if (i %% 500 == 0) message("  Processed ", i, "/", pb_total, " genes")
    }

    # jsonliteでシンプル&最小化 (auto_unbox = TRUE にすると無駄な配列カッコが消える)
    write_json(list(genes = names(expr_list), expr = expr_list), out_json_path, pretty = FALSE, auto_unbox = TRUE)
    message("  Saved: ", out_json_path, "  (", length(expr_list), " genes)")
    return(TRUE)
}

# ── 実行ロジック ────────────────────────────────────────

if (process_all) {
    # Individual samples (GSE145926)
    qs2_files <- list.files(OUTPUT_DIR, pattern = "^seurat_GSE145926_.*\\.qs2$", full.names = FALSE)
    qs2_files <- qs2_files[!grepl("merged", qs2_files)]
    for (f in qs2_files) {
        full_id <- sub("^seurat_", "", sub("\\.qs2$", "", f))
        process_file(file.path(OUTPUT_DIR, f), file.path(OUTPUT_DIR, paste0("expr_data_", full_id, ".json")), full_id)
    }
}

if (!is.null(sample_id)) {
    full_id <- paste0("GSE145926_", sample_id)
    process_file(
        file.path(OUTPUT_DIR, paste0("seurat_", full_id, ".qs2")),
        file.path(OUTPUT_DIR, paste0("expr_data_", full_id, ".json")), full_id
    )
}

if (!is.null(merged_id)) {
    # Merged study (lowercase usually)
    qs2_path <- file.path(OUTPUT_DIR, paste0("seurat_", merged_id, "_merged.qs2"))
    out_path <- file.path(OUTPUT_DIR, paste0("expr_data_", merged_id, "_merged.json"))
    process_file(qs2_path, out_path, paste0(merged_id, " (merged)"))
}

if (process_mvp) {
    qs2_path <- file.path(OUTPUT_DIR, "seurat_MVP_Sample.qs2")
    out_path <- file.path(OUTPUT_DIR, "expr_data_mvp.json") # app.js match
    process_file(qs2_path, out_path, "MVP_Sample")
}

if (!process_all && is.null(sample_id) && is.null(merged_id) && !process_mvp) {
    message("Usage:\n  Rscript generate_expr_data.R --sample C141\n  Rscript generate_expr_data.R --merged gse145926\n  Rscript generate_expr_data.R --mvp\n  Rscript generate_expr_data.R --all")
}
