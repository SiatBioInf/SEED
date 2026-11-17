#!/usr/bin/env Rscript

# Convert SEED h5ad exchange files to Seurat RDS files.
# Run this once before the downstream R/Seurat workflow.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(normalizePath(file.path(script_dir, "..", "01_config", "01_setup_paths.R"), mustWork = TRUE))
required_packages(c("Seurat"))

overwrite <- "--overwrite" %in% commandArgs(trailingOnly = TRUE)

convert_with_seuratdisk <- function(h5ad_path) {
  required_packages(c("SeuratDisk"))
  h5seurat_path <- sub("\\.h5ad$", ".h5seurat", h5ad_path, ignore.case = TRUE)
  SeuratDisk::Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
  SeuratDisk::LoadH5Seurat(h5seurat_path)
}

convert_with_zellkonverter <- function(h5ad_path) {
  required_packages(c("zellkonverter", "SingleCellExperiment", "SummarizedExperiment"))
  sce <- zellkonverter::readH5AD(h5ad_path)
  assay_names <- SummarizedExperiment::assayNames(sce)
  if (length(assay_names) == 0) {
    stop("No assay matrix found in h5ad: ", h5ad_path, call. = FALSE)
  }
  counts_assay <- if ("counts" %in% assay_names) "counts" else assay_names[1]
  data_assay <- if ("logcounts" %in% assay_names) "logcounts" else counts_assay
  Seurat::as.Seurat(sce, counts = counts_assay, data = data_assay)
}

convert_one <- function(h5ad_path, rds_path, label, overwrite = FALSE) {
  if (!file.exists(h5ad_path)) {
    message("Skip ", label, ": h5ad not found: ", h5ad_path)
    return(invisible(NULL))
  }
  if (file.exists(rds_path) && !overwrite) {
    message("Skip ", label, ": RDS already exists: ", rds_path)
    return(invisible(NULL))
  }

  message("Converting ", label, ": ", h5ad_path)
  obj <- tryCatch(
    convert_with_seuratdisk(h5ad_path),
    error = function(e) {
      message("SeuratDisk conversion failed for ", label, ": ", conditionMessage(e))
      message("Trying zellkonverter fallback for ", label, ".")
      convert_with_zellkonverter(h5ad_path)
    }
  )
  saveRDS(obj, rds_path)
  message("Saved ", label, " RDS: ", rds_path)
  invisible(rds_path)
}

convert_one(PATHS$sc_h5ad, PATHS$sc_rds, "SEED scRNA", overwrite)
convert_one(PATHS$bulk_h5ad, PATHS$bulk_rds, "SEED bulk RNA", overwrite)
convert_one(PATHS$pbmc_h5ad, PATHS$pbmc_rds, "non-CTC PBMC", overwrite)

message("Conversion finished.")
