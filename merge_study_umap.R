# ============================================================
#  merge_study_umap.R
#  Study全体のサンプルをMerge・Harmonyでバッチ補正し、
#  Study-level Merged UMAPデータをJSON出力する
#
#  Usage:
#    Rscript merge_study_umap.R --study gse145926
#
#  Run in: conda activate alew
#
#  入力: output/seurat_GSE145926_CXXX.qs2  (既存の個別処理済みファイル)
#  出力:
#    output/umap_data_gse145926_merged.json
#    output/top_clones_gse145926_merged.json
#    output/seurat_gse145926_merged.qs2
# ============================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(harmony) # Harmonyバッチ補正
    library(scRepertoire)
    library(tidyverse)
    library(jsonlite)
    library(qs2)
})

# ── コマンドライン引数 ────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
study_id <- "gse145926" # デフォルト
for (i in seq_along(args)) {
    if (args[i] == "--study" && i < length(args)) {
        study_id <- args[i + 1]
    }
}

message("======================================================")
message(" Merge Study UMAP")
message(" Study ID: ", study_id)
message("======================================================")

# ── Paths ─────────────────────────────────────────────────
BASE_DIR <- "/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database"
OUTPUT_DIR <- file.path(BASE_DIR, "output")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ── Study設定 ─────────────────────────────────────────────
STUDY_CONFIG <- list(
    gse145926 = list(
        samples = c(
            "C141", "C142", "C143", "C144", "C145",
            "C146", "C51", "C52", "C100", "C148",
            "C149", "C152"
        ),
        severity = c(
            C141 = "Moderate", C142 = "Moderate", C143 = "Moderate",
            C144 = "Severe", C145 = "Severe", C146 = "Severe",
            C51 = "Severe", C52 = "Moderate", C100 = "Moderate",
            C148 = "Severe", C149 = "Severe", C152 = "Severe"
        ),
        prefix = "GSE145926",
        full_name = "GSE145926 — Liao et al. Nature Medicine 2020"
    )
)

if (!study_id %in% names(STUDY_CONFIG)) {
    stop(
        "Unknown study_id: ", study_id,
        "\nAvailable: ", paste(names(STUDY_CONFIG), collapse = ", ")
    )
}

cfg <- STUDY_CONFIG[[study_id]]

# ── Step 1: 既存 qs2 から Seurat オブジェクトをロード ──
message("\n=== Step 1: Loading individual Seurat objects ===")

seurat_list <- list()
loaded_samples <- c()

for (sid in cfg$samples) {
    full_id <- paste0(cfg$prefix, "_", sid)
    qs2_path <- file.path(OUTPUT_DIR, paste0("seurat_", full_id, ".qs2"))

    if (!file.exists(qs2_path)) {
        message("  [SKIP] Not found: ", qs2_path)
        next
    }

    message("  Loading: ", qs2_path)
    seu <- tryCatch(
        qs_read(qs2_path),
        error = function(e) {
            message("  [ERROR] Failed to load: ", e$message)
            NULL
        }
    )
    if (is.null(seu)) next

    # サンプルメタデータを付与
    seu$sample_id <- sid
    seu$full_id <- full_id
    seu$severity <- cfg$severity[[sid]]
    seu$study_id <- study_id

    seurat_list[[sid]] <- seu
    loaded_samples <- c(loaded_samples, sid)
    message("    Cells: ", ncol(seu), "  Clusters: ", length(unique(seu$seurat_clusters)))
}

if (length(seurat_list) < 2) {
    stop(
        "Need at least 2 samples to merge. ",
        "Run process_geo_sample.R for each sample first."
    )
}

message(
    "\n  Loaded ", length(seurat_list), " samples: ",
    paste(loaded_samples, collapse = ", ")
)

# ── Step 2: Merge ─────────────────────────────────────────
message("\n=== Step 2: Merging Seurat objects ===")

# 最初のオブジェクト + 残りをmerge
seu_merged <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = loaded_samples, # バーコードに sample_id プレフィックス追加
    project = paste0(study_id, "_merged")
)

message("  Total cells after merge: ", ncol(seu_merged))
message(
    "  Samples in metadata: ",
    paste(sort(unique(seu_merged$sample_id)), collapse = ", ")
)

# ── Step 3: 標準前処理 (Mergeオブジェクト用) ──────────────
message("\n=== Step 3: Preprocessing merged object ===")

# Mergeした後は再度正規化 & 変数遺伝子選択が必要
seu_merged <- NormalizeData(seu_merged, verbose = FALSE)
seu_merged <- FindVariableFeatures(seu_merged, nfeatures = 3000, verbose = FALSE)
seu_merged <- ScaleData(seu_merged, verbose = FALSE)
seu_merged <- RunPCA(seu_merged, npcs = 30, verbose = FALSE)

message("  PCA done.")

# ── Step 4: Harmony バッチ補正 ────────────────────────────
message("\n=== Step 4: Harmony batch correction (by sample_id) ===")

# Harmonyがインストールされているか確認
if (!requireNamespace("harmony", quietly = TRUE)) {
    message("  [WARNING] harmony not installed. Falling back to no batch correction.")
    message("  Install with: conda run -n alew R -e \"install.packages('harmony')\"")
    # フォールバック: PCAをそのまま使ってUMAP
    seu_merged <- RunUMAP(seu_merged, dims = 1:30, verbose = FALSE)
    reduction_use <- "pca"
} else {
    seu_merged <- RunHarmony(
        seu_merged,
        group.by.vars = "sample_id", # sample_idでバッチ補正
        reduction = "pca",
        reduction.save = "harmony",
        verbose = FALSE
    )
    message("  Harmony done.")
    reduction_use <- "harmony"
}

# ── Step 5: UMAP & Clustering (Harmonyを使用) ─────────────
message("\n=== Step 5: UMAP & Clustering ===")

seu_merged <- RunUMAP(
    seu_merged,
    reduction = reduction_use,
    dims      = 1:30,
    verbose   = FALSE
)
seu_merged <- FindNeighbors(
    seu_merged,
    reduction = reduction_use,
    dims      = 1:30,
    verbose   = FALSE
)
# resolution 0.4 = 少し大きめのクラスタ (多サンプルMergeなので)
seu_merged <- FindClusters(seu_merged, resolution = 0.4, verbose = FALSE)

n_clusters <- length(unique(seu_merged$seurat_clusters))
message("  UMAP done. Clusters found: ", n_clusters)
message("  Cells per sample:")
print(table(seu_merged$sample_id))

# ── Step 6: Save Merged Seurat object ────────────────────
message("\n=== Step 6: Saving merged Seurat object ===")

out_qs2 <- file.path(OUTPUT_DIR, paste0("seurat_", study_id, "_merged.qs2"))
qs_save(seu_merged, out_qs2)
message("  Saved: ", out_qs2)

# ── Step 7: JSON出力 ──────────────────────────────────────
message("\n=== Step 7: Exporting JSON for WebApp ===")

umap_coords <- Embeddings(seu_merged, "umap") %>%
    as.data.frame() %>%
    setNames(c("umap_1", "umap_2"))

metadata <- seu_merged@meta.data %>%
    rownames_to_column("barcode") %>%
    select(
        barcode,
        seurat_clusters,
        sample_id, # ← サンプルフィルタリング用
        severity, # ← Severity色分け用
        any_of(c(
            "TCR_CTgene", "TCR_CTnt", "TCR_CTaa",
            "BCR_CTgene", "BCR_CTnt", "BCR_CTaa"
        ))
    ) %>%
    # seurat_clusters を文字列化
    mutate(seurat_clusters = as.character(seurat_clusters))

export_df <- bind_cols(metadata, umap_coords)

# JSON出力
out_json <- file.path(OUTPUT_DIR, paste0("umap_data_", study_id, "_merged.json"))
write_json(export_df, out_json, pretty = FALSE)
message("  Saved: ", out_json, "  (", nrow(export_df), " cells)")

# ── Step 8: Top Clones (Study全体 + Sample別) ────────────
message("\n=== Step 8: Exporting top clones ===")

get_top_clones_per_sample <- function(df, prefix, top_n = 20) {
    col_name <- paste0(prefix, "_CTgene")
    if (!col_name %in% names(df)) {
        return(NULL)
    }

    filtered <- df %>% filter(!is.na(!!sym(col_name)))
    if (nrow(filtered) == 0) {
        return(NULL)
    }

    # Study全体のTop clones
    overall <- filtered %>%
        count(!!sym(col_name), sort = TRUE) %>%
        head(top_n) %>%
        rename(clonotype_id = !!sym(col_name), count = n) %>%
        mutate(
            clonotype_id = as.character(clonotype_id),
            count        = as.integer(count),
            type         = prefix,
            sample_id    = "ALL"
        )

    # サンプル別Top clones
    per_sample <- filtered %>%
        group_by(sample_id) %>%
        count(!!sym(col_name), sort = TRUE) %>%
        slice_head(n = 10) %>%
        ungroup() %>%
        rename(clonotype_id = !!sym(col_name), count = n) %>%
        mutate(
            clonotype_id = as.character(clonotype_id),
            count        = as.integer(count),
            type         = prefix
        )

    bind_rows(overall, per_sample)
}

tcr_clones <- get_top_clones_per_sample(export_df, "TCR", top_n = 20)
bcr_clones <- get_top_clones_per_sample(export_df, "BCR", top_n = 20)

clone_list <- Filter(Negate(is.null), list(tcr_clones, bcr_clones))
if (length(clone_list) > 0) {
    top_clones <- bind_rows(clone_list)
} else {
    top_clones <- data.frame(
        clonotype_id = character(0),
        count        = integer(0),
        type         = character(0),
        sample_id    = character(0)
    )
}

out_clones <- file.path(OUTPUT_DIR, paste0("top_clones_", study_id, "_merged.json"))
write_json(top_clones, out_clones, pretty = TRUE)
message("  Saved: ", out_clones)

# ── Step 9: Sample summary (Webアプリ用) ─────────────────
message("\n=== Step 9: Generating sample summary ===")

# TCR列が存在する場合のみカウント
if ("TCR_CTgene" %in% names(export_df)) {
    sample_summary <- export_df %>%
        group_by(sample_id, severity) %>%
        summarise(
            n_cells    = n(),
            n_clusters = n_distinct(seurat_clusters),
            n_tcr      = sum(!is.na(TCR_CTgene)),
            .groups    = "drop"
        ) %>%
        arrange(severity, sample_id)
} else {
    sample_summary <- export_df %>%
        group_by(sample_id, severity) %>%
        summarise(
            n_cells    = n(),
            n_clusters = n_distinct(seurat_clusters),
            n_tcr      = 0L,
            .groups    = "drop"
        ) %>%
        arrange(severity, sample_id)
}

out_summary <- file.path(OUTPUT_DIR, paste0("sample_summary_", study_id, "_merged.json"))
write_json(sample_summary, out_summary, pretty = TRUE)
message("  Saved: ", out_summary)
print(sample_summary)

# ── 完了メッセージ ────────────────────────────────────────
message("\n======================================================")
message(" DONE: ", study_id, " Merged UMAP")
message("  Total cells  : ", nrow(export_df))
message("  Clusters     : ", n_clusters)
message("  Samples      : ", paste(loaded_samples, collapse = ", "))
message("  Output files :")
message("    ", out_json)
message("    ", out_clones)
message("    ", out_summary)
message("    ", out_qs2)
message("======================================================")
message("\nNext step: Update studies.json to add merged viewer entry")
message("  Then open: http://<server>:8080/webapp/viewer.html?study=gse145926")
