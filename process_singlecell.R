# Load required libraries
library(Seurat)
library(scRepertoire)
library(tidyverse)
library(jsonlite)
library(qs2)

# Set working directory and paths
data_dir <- "/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database/data"
output_dir <- "/home/masakazu/Data/260212_alew_singlecell_tcrbcr_database/output"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define file paths
h5_file <- file.path(data_dir, "10k_5p_Human_diseased_PBMC_ALL_Fresh_count_filtered_feature_bc_matrix.h5")
tcr_file <- file.path(data_dir, "10k_5p_Human_diseased_PBMC_ALL_Fresh_vdj_t_filtered_contig_annotations.csv")
bcr_file <- file.path(data_dir, "10k_5p_Human_diseased_PBMC_ALL_Fresh_vdj_b_filtered_contig_annotations.csv")

# --- Step 1: Process Gene Expression Data ---
message("Processing Gene Expression Data...")
seu.data <- Read10X_h5(h5_file)
seu <- CreateSeuratObject(counts = seu.data, project = "MVP_Sample")

# Standard Seurat pipeline
# Rename cells to match scRepertoire format (Sample_ID_Barcode)
# combined_tcr/bcr will have prefixes "MVP_Sample_" containing the sample and ID
seu <- RenameCells(seu, new.names = paste0("MVP_Sample_", colnames(seu)))

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)

# --- Step 2: Process VDJ Data (TCR/BCR) ---
message("Processing VDJ Data...")
tcr <- read.csv(tcr_file)
bcr <- read.csv(bcr_file)

# Combine TCR and BCR using scRepertoire
# We treat them as a list for combineTCR/combineBCR
combined_tcr <- combineTCR(list(tcr), samples = "MVP", ID = "Sample")
combined_bcr <- combineBCR(list(bcr), samples = "MVP", ID = "Sample")

# Combine Expression with TCR/BCR
# Note: scRepertoire's combineExpression adds clonotype info to the Seurat object metadata
seu <- combineExpression(combined_tcr, seu, cloneCall = "gene", proportion = FALSE, cloneSize = "Gene+Nt", addLabel = FALSE)
# Combine BCR as well? combineExpression is designed for one, checking docs...
# Usually we run it once. If we want both, we might need to merge or run sequentially.
# For this MVP, let's focus on TCR first or see if we can add BCR columns manually or run it again.
# Running combineExpression again might overwrite or add new columns. Let's try adding BCR as well.
# It seems `combineExpression` attaches to `@meta.data`. Let's try to add BCR info carefully.
# ACTUALLY, a better approach for MVP to show *clonal* relationships might be just one (TCR or BCR) or merging the contigs list?
# `combineTCR` and `combineBCR` return lists of data frames.
# Let's try to merge the contig lists if possible, or just focus on adding them one by one.
# For simplicity in MVP, let's add TCR first, then BCR, but renaming columns to avoid overwrite if needed.
# However, `combineExpression` is robust. Let's try adding BCR second.
# Warning: The second call might overwrite 'CTgene', 'CTnt', etc.
# Let's inspect the column names after first call? We can't interactively.
# Strategy: Run for TCR. Rename columns (prefix 'TCR_'). Run for BCR. Rename columns (prefix 'BCR_').

# Re-loading to be safe for a clean merge strategy
# Actually, scRepertoire 2.0+ might handle this better, but let's stick to standard flow.
# Let's just create a combined list of contigs? No, they are different structure (TRA/TRB vs IGH/IGL).
# Let's add TCR info.
seu <- combineExpression(combined_tcr, seu, cloneCall = "gene")
# Rename columns to protect them
tcr_cols <- c("CTgene", "CTnt", "CTaa", "CTstrict", "cloneType")
for (col in tcr_cols) {
  if (col %in% colnames(seu@meta.data)) {
    seu@meta.data[[paste0("TCR_", col)]] <- seu@meta.data[[col]]
    seu@meta.data[[col]] <- NULL # Remove original to avoid confusion
  }
}

# Add BCR info
seu <- combineExpression(combined_bcr, seu, cloneCall = "gene")
# Rename columns
bcr_cols <- c("CTgene", "CTnt", "CTaa", "CTstrict", "cloneType")
for (col in bcr_cols) {
  if (col %in% colnames(seu@meta.data)) {
    seu@meta.data[[paste0("BCR_", col)]] <- seu@meta.data[[col]]
    seu@meta.data[[col]] <- NULL
  }
}


# --- Step 3: Export Lightweight Data for Web ---
message("Exporting Data for Web...")

# qs2形式でSeuratオブジェクト全体を保存
out_qs2 <- file.path(output_dir, "seurat_MVP_Sample.qs2")
qs_save(seu, out_qs2)
message("Saved (qs2): ", out_qs2)

# 1. UMAP Coordinates + Metadata
umap_coords <- Embeddings(seu, "umap") %>% as.data.frame()
metadata <- seu@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode, seurat_clusters, starts_with("TCR_"), starts_with("BCR_"))

export_df <- bind_cols(metadata, umap_coords)

# Save as JSON
write_json(export_df, file.path(output_dir, "umap_data.json"), pretty = FALSE)

# 2. Clonotype Summary (Top 10 TCR and BCR)
# Calculate frequencies
get_top_clones <- function(df, prefix) {
  col_name <- paste0(prefix, "_CTgene")
  if (!col_name %in% names(df)) {
    return(NULL)
  }

  df %>%
    filter(!is.na(!!sym(col_name))) %>%
    count(!!sym(col_name), sort = TRUE) %>%
    head(10) %>%
    rename(clonotype_id = !!sym(col_name), count = n) %>%
    mutate(type = prefix)
}

top_tcr <- get_top_clones(export_df, "TCR")
top_bcr <- get_top_clones(export_df, "BCR")

top_clones <- bind_rows(top_tcr, top_bcr)
write_json(top_clones, file.path(output_dir, "top_clones.json"), pretty = TRUE)

# 3. Export expr_data JSON for Gene Feature Plot
message("Exporting expr_data JSON...")
N_GENES <- 3000
var_features <- VariableFeatures(seu)
if (length(var_features) == 0) {
  seu <- FindVariableFeatures(seu, nfeatures = N_GENES)
  var_features <- VariableFeatures(seu)
}
top_genes <- head(var_features, N_GENES)
expr_mat <- GetAssayData(seu, layer = "data")[top_genes, , drop = FALSE]
cell_barcodes <- export_df$barcode

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

out_expr <- file.path(output_dir, "expr_data_mvp.json")
write_json(list(genes = names(expr_list), expr = expr_list), out_expr, pretty = FALSE)
message("Saved: ", out_expr)

message("Done! Output files saved to: ", output_dir)
message()
message("To reload the Seurat object in R:")
message("  library(qs2)")
message("  seu <- qs_read(file.path(output_dir, 'seurat_MVP_Sample.qs2'))")
