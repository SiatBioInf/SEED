#!/usr/bin/env Rscript

# Calculate 31 SEED state scores and median-based High/Low group labels.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(normalizePath(file.path(script_dir, "..", "01_config", "01_setup_paths.R"), mustWork = TRUE))
required_packages(c("Seurat", "dplyr", "Matrix"))

input_rds <- file.path(PATHS$output, "02_SEED_scRNA_preprocessed.rds")
if (!file.exists(input_rds)) {
  input_rds <- PATHS$sc_rds
}
output_rds <- file.path(PATHS$output, "03_SEED_scRNA_scored.rds")

load_gene_sets_from_csv <- function(csv_path, expected_scores = STATE_THEMES$score) {
  if (!file.exists(csv_path)) {
    stop("Curated gene-set CSV not found: ", csv_path, call. = FALSE)
  }

  gene_set_table <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("State_score", "Gene_symbols")
  missing_cols <- setdiff(required_cols, colnames(gene_set_table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in curated gene-set CSV: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  gene_set_table <- gene_set_table[
    !is.na(gene_set_table[["State_score"]]) & gene_set_table[["State_score"]] != "",
    ,
    drop = FALSE
  ]
  gene_set_table <- gene_set_table[!duplicated(gene_set_table[["State_score"]]), , drop = FALSE]

  gene_sets <- lapply(gene_set_table[["Gene_symbols"]], function(x) {
    clean_genes(unlist(strsplit(x, ";", fixed = TRUE)))
  })
  names(gene_sets) <- gene_set_table[["State_score"]]

  missing_scores <- setdiff(expected_scores, names(gene_sets))
  if (length(missing_scores) > 0) {
    stop("Curated gene-set CSV is missing state scores: ", paste(missing_scores, collapse = ", "), call. = FALSE)
  }

  gene_sets <- gene_sets[expected_scores]
  empty_sets <- names(gene_sets)[vapply(gene_sets, length, integer(1)) < 3]
  if (length(empty_sets) > 0) {
    stop("Too few genes in curated gene sets: ", paste(empty_sets, collapse = ", "), call. = FALSE)
  }

  long_table <- do.call(
    rbind,
    lapply(names(gene_sets), function(score_name) {
      rows <- gene_set_table[gene_set_table[["State_score"]] == score_name, , drop = FALSE]
      data.frame(
        theme = STATE_THEMES$theme[match(score_name, STATE_THEMES$score)],
        state_score = score_name,
        state_name = rows[["State_name"]][1],
        gene_set_source = rows[["Gene_set_source"]][1],
        gene = gene_sets[[score_name]],
        stringsAsFactors = FALSE
      )
    })
  )
  write.csv(long_table, file.path(PATHS$release_output, "SEED_curated_gene_sets_long.csv"), row.names = FALSE)

  gene_sets
}

obj <- require_rds(input_rds, "preprocessed or raw SEED scRNA object")
gene_sets <- load_gene_sets_from_csv(PATHS$gene_sets_csv)
obj <- set_default_assay(obj)

score_one_set <- function(object, score_name, genes) {
  genes_use <- intersect(clean_genes(genes), rownames(object))
  if (length(genes_use) < 3) {
    warning("Too few genes found for ", score_name, ": ", length(genes_use))
    object[[score_name]] <- NA_real_
    return(object)
  }
  tmp_name <- paste0(score_name, "_tmp")
  object <- AddModuleScore(object, features = list(genes_use), name = tmp_name, search = FALSE)
  generated_col <- paste0(tmp_name, "1")
  object[[score_name]] <- object@meta.data[[generated_col]]
  object@meta.data[[generated_col]] <- NULL
  object
}

for (score_name in names(gene_sets)) {
  message("Scoring ", score_name)
  obj <- score_one_set(obj, score_name, gene_sets[[score_name]])
  med <- median(obj@meta.data[[score_name]], na.rm = TRUE)
  group_col <- paste0(score_name, "_Group")
  obj[[group_col]] <- factor(ifelse(obj@meta.data[[score_name]] > med, "High", "Low"), levels = c("Low", "High"))
}

saveRDS(obj, output_rds)
write.csv(obj@meta.data, file.path(PATHS$release_output, "scRNA_metadata_with_state_scores.csv"))
message("Saved: ", output_rds)
