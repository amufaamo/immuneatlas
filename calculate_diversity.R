# ============================================================
#  calculate_diversity.R
#  各サンプルのクローン多様性指標（Shannon, Simpson 等）を計算する
#  Usage:
#    conda activate alew
#    Rscript calculate_diversity.R
# ============================================================

library(Seurat)
library(scRepertoire)
library(tidyverse)
library(jsonlite)
library(qs2)

BASE_DIR <- "/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database"
OUTPUT_DIR <- file.path(BASE_DIR, "output")

# output 取得
qs2_files <- list.files(OUTPUT_DIR, pattern = "^seurat_GSE145926_.*\\.qs2$", full.names = TRUE)
# 個別サンプルのみ（mergedは別階層または別名にする）
qs2_files <- qs2_files[!grepl("merged", qs2_files)]

results <- list()

for (f in qs2_files) {
    full_id <- sub("^seurat_", "", sub("\\.qs2$", "", basename(f)))
    message("Processing: ", full_id)

    seu <- qs_read(f)

    # TCR多様性
    tcr_df <- seu@meta.data %>% filter(!is.na(TCR_CTgene))

    if (nrow(tcr_df) > 10) {
        # Clonal diversity using scRepertoire (internal helper or manual)
        # scRepertoire::clonotypeDiversity は Seurat オブジェクトを直接受け取れるが、
        # ここではシンプルに計算する

        counts <- table(tcr_df$TCR_CTgene)
        props <- counts / sum(counts)

        shannon <- -sum(props * log(props))
        simpson <- sum(props^2) # Higher = less diverse, usually 1-D is used
        inv_simpson <- 1 / simpson

        div <- list(
            sample_id = full_id,
            n_cells = ncol(seu),
            n_tcr = nrow(tcr_df),
            shannon = round(shannon, 4),
            simpson = round(simpson, 4),
            inv_simpson = round(inv_simpson, 4)
        )
    } else {
        div <- list(
            sample_id = full_id,
            n_cells = ncol(seu),
            n_tcr = nrow(tcr_df),
            shannon = NA,
            simpson = NA,
            inv_simpson = NA
        )
    }
    results[[full_id]] <- div
}

# MVPサンプル
mvp_f <- file.path(OUTPUT_DIR, "seurat_MVP_Sample.qs2")
if (file.exists(mvp_f)) {
    message("Processing: MVP_Sample")
    seu <- qs_read(mvp_f)
    tcr_df <- seu@meta.data %>% filter(!is.na(TCR_CTgene))
    if (nrow(tcr_df) > 10) {
        counts <- table(tcr_df$TCR_CTgene)
        props <- counts / sum(counts)
        shannon <- -sum(props * log(props))
        simpson <- sum(props^2)
        inv_simpson <- 1 / simpson
        results[["MVP_Sample"]] <- list(
            sample_id = "MVP_Sample",
            n_cells = ncol(seu),
            n_tcr = nrow(tcr_df),
            shannon = round(shannon, 4),
            simpson = round(simpson, 4),
            inv_simpson = round(inv_simpson, 4)
        )
    }
}

# Integrated / Merged Study
merged_f <- file.path(OUTPUT_DIR, "seurat_gse145926_merged.qs2")
if (file.exists(merged_f)) {
    message("Processing: gse145926_merged")
    seu <- qs_read(merged_f)
    # merged では TCR_CTgene などのカラム名が JoinLayers 後に残っているか確認が必要だが、
    # 通常はそのまま残る。
    tcr_df <- seu@meta.data %>% filter(!is.na(TCR_CTgene))

    if (nrow(tcr_df) > 10) {
        counts <- table(tcr_df$TCR_CTgene)
        props <- counts / sum(counts)
        shannon <- -sum(props * log(props))
        simpson <- sum(props^2)
        inv_simpson <- 1 / simpson
        results[["gse145926_merged"]] <- list(
            sample_id = "gse145926_merged",
            n_cells = ncol(seu),
            n_tcr = nrow(tcr_df),
            shannon = round(shannon, 4),
            simpson = round(simpson, 4),
            inv_simpson = round(inv_simpson, 4)
        )
    } else {
        results[["gse145926_merged"]] <- list(
            sample_id = "gse145926_merged",
            n_cells = ncol(seu),
            n_tcr = nrow(tcr_df),
            shannon = NA, simpson = NA, inv_simpson = NA
        )
    }
}

write_json(results, file.path(OUTPUT_DIR, "diversity_summary.json"), auto_unbox = TRUE, pretty = TRUE)
message("Saved: diversity_summary.json")
