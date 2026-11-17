#!/usr/bin/env Rscript

# Merge harmonized metadata and run a standard scRNA-seq preprocessing workflow.
# Multi-resolution cluster selection and ShinyCell generation are intentionally omitted.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(normalizePath(file.path(script_dir, "..", "01_config", "01_setup_paths.R"), mustWork = TRUE))
required_packages(c("Seurat", "dplyr", "readxl", "ggplot2", "patchwork"))

output_rds <- file.path(PATHS$output, "02_SEED_scRNA_preprocessed.rds")

# Optional external metadata table. The released SEED h5ad/RDS already stores
# harmonized metadata, so this is normally left as NULL.
extra_metadata_xlsx <- NULL

obj <- require_rds(PATHS$sc_rds, "SEED scRNA")

if (!is.null(extra_metadata_xlsx) && file.exists(extra_metadata_xlsx)) {
  meta_use <- readxl::read_excel(extra_metadata_xlsx) |> as.data.frame()
  if (!"CTC_ID" %in% colnames(meta_use)) {
    stop("Metadata table must contain a CTC_ID column.", call. = FALSE)
  }
  rownames(meta_use) <- meta_use$CTC_ID
  meta_use$CTC_ID <- NULL
  common_cells <- intersect(colnames(obj), rownames(meta_use))
  obj <- subset(obj, cells = common_cells)
  obj@meta.data <- cbind(obj@meta.data, meta_use[colnames(obj), , drop = FALSE])
}

obj <- set_default_assay(obj)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)

if ("dataset" %in% colnames(obj@meta.data) && requireNamespace("harmony", quietly = TRUE)) {
  obj <- harmony::RunHarmony(obj, group.by.vars = "dataset", dims.use = 1:30, verbose = TRUE)
  obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony")
  obj <- FindNeighbors(obj, dims = 1:30, reduction = "harmony")
} else {
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca")
  obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
}

if ("cancer" %in% colnames(obj@meta.data)) {
  obj$cancer <- factor(obj$cancer, levels = names(CANCER_COLORS))
}

saveRDS(obj, output_rds)
write.csv(obj@meta.data, file.path(PATHS$release_output, "scRNA_metadata_after_preprocessing.csv"))
message("Saved: ", output_rds)
