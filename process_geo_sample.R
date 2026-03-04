# ============================================================
#  process_geo_sample.R  (汎用版)
#  GSE145926 各サンプルの処理スクリプト
#  Usage:
#    Rscript process_geo_sample.R --sample C141
#    Rscript process_geo_sample.R --sample C142
#  Run in: conda activate alew
# ============================================================

library(Seurat)
library(scRepertoire)
library(tidyverse)
library(jsonlite)
library(qs2)

# ── コマンドライン引数パース ────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

# --sample XXXX を探す
sample_id <- NULL
for (i in seq_along(args)) {
    if (args[i] == "--sample" && i < length(args)) {
        sample_id <- args[i + 1]
    }
}

# 引数がない場合はデフォルト (後方互換)
if (is.null(sample_id)) {
    sample_id <- "C141"
    message("No --sample argument provided. Using default: ", sample_id)
}

# ── Paths ───────────────────────────────────────────────────
BASE_DIR <- "/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database"
DATA_DIR <- file.path(BASE_DIR, "data", paste0("GSE145926_", sample_id))
OUTPUT_DIR <- file.path(BASE_DIR, "output")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

message("======================================================")
message(" Processing sample: ", sample_id)
message(" Data dir        : ", DATA_DIR)
message("======================================================")

if (!dir.exists(DATA_DIR)) {
    stop(
        "Data directory not found: ", DATA_DIR,
        "\nRun download_gse145926.sh first."
    )
}

h5_file <- file.path(DATA_DIR, "filtered_feature_bc_matrix.h5")
tcr_file <- file.path(DATA_DIR, "filtered_contig_annotations_tcr.csv")

if (!file.exists(h5_file)) {
    stop("H5 file not found: ", h5_file)
}

# ── Step 1: GEX ─────────────────────────────────────────────
message("=== Step 1: Loading GEX ===")

SEU_PREFIX <- paste0("GEO_GSE145926_", sample_id)
FULL_ID <- paste0("GSE145926_", sample_id)

seu.data <- Read10X_h5(h5_file)
if (is.list(seu.data)) {
    message("Multi-modal H5 detected, using 'Gene Expression' slot")
    seu.data <- seu.data[["Gene Expression"]]
}

seu <- CreateSeuratObject(
    counts = seu.data, project = FULL_ID,
    min.cells = 3, min.features = 200
)
seu <- RenameCells(seu, new.names = paste0(SEU_PREFIX, "_", colnames(seu)))

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30)
seu <- RunUMAP(seu, dims = 1:20)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.5)

message("Cells after QC: ", ncol(seu))

# ── Step 2: VDJ ─────────────────────────────────────────────
message("=== Step 2: Loading VDJ ===")

if (!is.null(tcr_file) && file.exists(tcr_file)) {
    tcr_raw <- read.csv(tcr_file)
    combined_tcr <- combineTCR(list(tcr_raw), samples = "GEO", ID = FULL_ID)
    seu <- combineExpression(combined_tcr, seu, cloneCall = "gene")
    tcr_cols <- c("CTgene", "CTnt", "CTaa", "CTstrict", "cloneType")
    for (col in tcr_cols) {
        if (col %in% colnames(seu@meta.data)) {
            seu@meta.data[[paste0("TCR_", col)]] <- seu@meta.data[[col]]
            seu@meta.data[[col]] <- NULL
        }
    }
    message("TCR combined. TCR+ cells: ", sum(!is.na(seu@meta.data$TCR_CTgene)))
} else {
    message("No TCR file found — filling with NA.")
    for (col in c("TCR_CTgene", "TCR_CTnt", "TCR_CTaa")) {
        seu@meta.data[[col]] <- NA
    }
}

# BCR — このデータセットはなし
for (col in c("BCR_CTgene", "BCR_CTnt", "BCR_CTaa")) {
    seu@meta.data[[col]] <- NA
}

# ── Step 3: Save Seurat as qs2 ──────────────────────────────
message("=== Step 3: Saving Seurat object as qs2 ===")
out_qs2 <- file.path(OUTPUT_DIR, paste0("seurat_", FULL_ID, ".qs2"))
qs_save(seu, out_qs2)
message("Saved (qs2): ", out_qs2)

# ── Step 4: Export JSON for Web ─────────────────────────────
message("=== Step 4: Exporting JSON ===")

umap_coords <- Embeddings(seu, "umap") %>% as.data.frame()
metadata <- seu@meta.data %>%
    rownames_to_column("barcode") %>%
    select(
        barcode, seurat_clusters,
        any_of(c(
            "TCR_CTgene", "TCR_CTnt", "TCR_CTaa",
            "BCR_CTgene", "BCR_CTnt", "BCR_CTaa"
        ))
    )

export_df <- bind_cols(metadata, umap_coords)

out_json <- file.path(OUTPUT_DIR, paste0("umap_data_", FULL_ID, ".json"))
out_clones <- file.path(OUTPUT_DIR, paste0("top_clones_", FULL_ID, ".json"))

write_json(export_df, out_json, pretty = FALSE)
message("Saved: ", out_json, "  (", nrow(export_df), " cells)")

# Top clones
get_top_clones <- function(df, prefix) {
    col_name <- paste0(prefix, "_CTgene")
    if (!col_name %in% names(df)) {
        return(NULL)
    }
    filtered <- df %>% filter(!is.na(!!sym(col_name)))
    if (nrow(filtered) == 0) {
        return(NULL)
    }
    filtered %>%
        count(!!sym(col_name), sort = TRUE) %>%
        head(10) %>%
        rename(clonotype_id = !!sym(col_name), count = n) %>%
        mutate(
            clonotype_id = as.character(clonotype_id),
            count = as.integer(count),
            type = prefix
        )
}

tcr_clones <- get_top_clones(export_df, "TCR")
bcr_clones <- get_top_clones(export_df, "BCR")

# NULLを除外してからbind_rows (型不一致エラー回避)
clone_list <- Filter(Negate(is.null), list(tcr_clones, bcr_clones))
if (length(clone_list) > 0) {
    top_clones <- bind_rows(clone_list)
} else {
    top_clones <- data.frame(
        clonotype_id = character(0),
        count = integer(0),
        type = character(0)
    )
}
write_json(top_clones, out_clones, pretty = TRUE)
message("Saved: ", out_clones)

# ── Step 5: Export expr_data JSON for Gene Feature Plot ─────
message("=== Step 5: Exporting expr_data JSON (Gene Feature Plot) ===")

# 高分散遺伝子 上位 N 個に絞る（ファイルサイズ・速度のため）
N_GENES <- 3000
var_features <- VariableFeatures(seu)
if (length(var_features) == 0) {
    seu <- FindVariableFeatures(seu, nfeatures = N_GENES)
    var_features <- VariableFeatures(seu)
}
top_genes <- head(var_features, N_GENES)

# log-normalized マトリクス (data slot) を取得
expr_mat <- GetAssayData(seu, layer = "data")[top_genes, , drop = FALSE]

# セルのバーコードを export_df のbarcode と対応させる
# export_df$barcode は "PREFIX_BARCODE-1" 形式
cell_barcodes <- export_df$barcode

# スパース辞書形式に変換: 0 は省略
# expr_list の構造: list(gene1 = list(barcode1 = val, ...), ...)
expr_list <- list()
pb_total <- length(top_genes)
for (i in seq_along(top_genes)) {
    gene <- top_genes[i]
    row_data <- expr_mat[i, ]
    nonzero_idx <- which(row_data > 0)
    if (length(nonzero_idx) == 0) next
    vals_nonzero <- as.numeric(row_data[nonzero_idx])
    vals_nonzero <- round(vals_nonzero, 3)
    names(vals_nonzero) <- cell_barcodes[nonzero_idx]
    expr_list[[gene]] <- as.list(vals_nonzero)
    if (i %% 500 == 0) message("  Processed ", i, "/", pb_total, " genes")
}

out_expr <- file.path(OUTPUT_DIR, paste0("expr_data_", FULL_ID, ".json"))
write_json(list(genes = names(expr_list), expr = expr_list), out_expr, pretty = FALSE)
message("Saved: ", out_expr, "  (", length(expr_list), " genes with non-zero expression)")

message("=== Done: ", sample_id, " ===")
message("  seurat_qs2 : ", out_qs2)
message("  umap_data  : ", out_json)
message("  top_clones : ", out_clones)
message("  expr_data  : ", out_expr)
