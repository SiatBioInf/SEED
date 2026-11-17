#!/usr/bin/env Rscript

# Inspect the curated gene sets without downloading the full atlas.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0) stop("Run this file with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", file_arg[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
gene_set_file <- file.path(
  repo_root, "data",
  "Supplementary Table S3. Curated gene sets and source references for the 31 functional states.csv"
)

sets <- read.csv(gene_set_file, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c("Theme", "State_score", "Gene_symbols", "Number_of_genes")
missing <- setdiff(needed, names(sets))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

sets <- sets[!is.na(sets$State_score) & nzchar(sets$State_score), , drop = FALSE]
if (anyDuplicated(sets$State_score)) stop("Duplicate state names in gene-set table.", call. = FALSE)
gene_counts <- vapply(strsplit(sets$Gene_symbols, ";", fixed = TRUE), length, integer(1))
if (!identical(as.integer(sets$Number_of_genes), gene_counts)) {
  stop("Declared gene counts do not match the gene symbols.", call. = FALSE)
}

cat(nrow(sets), "functional states in", length(unique(sets$Theme)), "themes\n")
cat("Genes per state: median", median(gene_counts), "; range", min(gene_counts), "to", max(gene_counts), "\n")
print(as.data.frame(table(Theme = sets$Theme)), row.names = FALSE)
