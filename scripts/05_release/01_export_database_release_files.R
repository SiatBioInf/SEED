#!/usr/bin/env Rscript

# Export SEED database release tables and file manifests.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(normalizePath(file.path(script_dir, "..", "01_config", "01_setup_paths.R"), mustWork = TRUE))
required_packages(c("Seurat", "dplyr", "fs", "writexl"))

release_dir <- PATHS$release_output
scored_rds <- file.path(PATHS$output, "03_SEED_scRNA_scored.rds")
if (!file.exists(scored_rds)) {
  scored_rds <- file.path(PATHS$output, "04_SEED_scRNA_scored.rds")
}

if (file.exists(scored_rds)) {
  sc <- readRDS(scored_rds)
  meta <- sc@meta.data
  write.csv(meta, file.path(release_dir, "SEED_scRNA_metadata.csv"))

  score_cols <- intersect(STATE_THEMES$score, colnames(meta))
  group_cols <- paste0(score_cols, "_Group")
  group_cols <- group_cols[group_cols %in% colnames(meta)]
  write.csv(meta[, score_cols, drop = FALSE], file.path(release_dir, "SEED_scRNA_state_scores.csv"))
  write.csv(meta[, group_cols, drop = FALSE], file.path(release_dir, "SEED_scRNA_state_groups.csv"))

  if ("cancer" %in% colnames(meta)) {
    write.csv(dplyr::count(meta, cancer, name = "cell_n"), file.path(release_dir, "SEED_scRNA_cancer_cell_count.csv"), row.names = FALSE)
  }
}

bulk_rds <- PATHS$bulk_rds
if (file.exists(bulk_rds)) {
  bulk <- readRDS(bulk_rds)
  write.csv(bulk@meta.data, file.path(release_dir, "SEED_bulkRNA_metadata.csv"))
}

make_manifest <- function(path, label) {
  if (!dir.exists(path)) return(data.frame())
  files <- fs::dir_ls(path = path, recurse = TRUE, type = "file")
  if (length(files) == 0) return(data.frame())
  data.frame(
    source = label,
    file_name = fs::path_file(files),
    relative_path = fs::path_rel(files, start = PATHS$root),
    extension = fs::path_ext(files),
    size_kb = round(as.numeric(fs::file_info(files)$size) / 1024, 2),
    modified_time = as.character(fs::file_info(files)$modification_time),
    stringsAsFactors = FALSE
  )
}

manifest <- dplyr::bind_rows(
  make_manifest(PATHS$confidence_output, "confidence_validation"),
  make_manifest(PATHS$sc_state_output, "scRNA_state_analysis"),
  make_manifest(PATHS$sc_cancer_output, "scRNA_cancer_vs_others"),
  make_manifest(PATHS$bulk_cancer_output, "bulkRNA_cancer_vs_others"),
  make_manifest(file.path(PATHS$root, "score_analysis_output-done"), "existing_state_results"),
  make_manifest(file.path(PATHS$root, "DEG"), "existing_cancer_deg_results"),
  make_manifest(file.path(PATHS$root, "browse_image"), "existing_browse_images")
)

write.csv(manifest, file.path(release_dir, "SEED_database_result_file_manifest.csv"), row.names = FALSE)
writexl::write_xlsx(list(file_manifest = manifest), path = file.path(release_dir, "SEED_database_result_file_manifest.xlsx"))

message("Release tables exported to: ", release_dir)
