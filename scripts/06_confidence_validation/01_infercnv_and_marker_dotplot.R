#!/usr/bin/env Rscript

# Biological confidence assessment for SEED CTCs.
# Outputs:
# 1) expression-derived inferCNV-like heatmap and CNV proxy score;
# 2) marker-gene dotplot for Experimental CTC, Computational CTC and PBMC.
#
# Official infercnv package execution is provided as an optional block at the end
# and is disabled by default because it is slower and needs more memory.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(normalizePath(file.path(script_dir, "..", "01_config", "01_setup_paths.R"), mustWork = TRUE))
required_packages(c("Seurat", "Matrix", "dplyr", "tidyr", "ggplot2", "scales"))

out_dir <- PATHS$confidence_output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ctc_rds <- file.path(PATHS$output, "03_SEED_scRNA_scored.rds")
if (!file.exists(ctc_rds)) {
  ctc_rds <- file.path(PATHS$output, "04_SEED_scRNA_scored.rds")
}
if (!file.exists(ctc_rds)) {
  ctc_rds <- PATHS$sc_rds
}

ctc <- require_rds(ctc_rds, "SEED scRNA")
pbmc <- require_rds(PATHS$pbmc_rds, "non-CTC PBMC")

ctc <- set_default_assay(ctc)
pbmc <- set_default_assay(pbmc)

get_count_matrix <- function(obj) {
  assay_use <- DefaultAssay(obj)
  mat <- tryCatch(
    GetAssayData(obj, assay = assay_use, layer = "counts"),
    error = function(e) GetAssayData(obj, assay = assay_use, slot = "counts")
  )
  if (!inherits(mat, "sparseMatrix")) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
  }
  mat
}

log_normalize_counts <- function(mat) {
  lib <- Matrix::colSums(mat)
  lib[lib <= 0] <- 1
  norm <- mat %*% Matrix::Diagonal(x = 1 / lib)
  norm <- norm * 1e4
  norm@x <- log1p(norm@x)
  norm
}

parse_gencode_gene_order <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    stop("Gene annotation GTF not found: ", gtf_path, call. = FALSE)
  }
  chrom_order <- c(as.character(1:22), "X", "Y")
  chrom_rank <- stats::setNames(seq_along(chrom_order), chrom_order)

  con <- gzfile(gtf_path, open = "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  parts <- parts[vapply(parts, length, integer(1)) >= 9]
  parts <- parts[vapply(parts, function(x) x[[3]] == "gene", logical(1))]

  rows <- lapply(parts, function(x) {
    chrom <- sub("^chr", "", x[[1]])
    if (!chrom %in% chrom_order) return(NULL)
    attrs <- x[[9]]
    gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", attrs)
    gene_type <- if (grepl('gene_type "', attrs)) sub('.*gene_type "([^"]+)".*', "\\1", attrs) else ""
    if (identical(gene_name, attrs) || is.na(gene_name) || gene_name == "") return(NULL)
    data.frame(
      gene = gene_name,
      chrom = chrom,
      start = as.integer(x[[4]]),
      end = as.integer(x[[5]]),
      gene_type = gene_type,
      stringsAsFactors = FALSE
    )
  })
  genes <- dplyr::bind_rows(rows)
  genes$priority <- dplyr::case_when(
    genes$gene_type == "protein_coding" ~ 0,
    genes$gene_type == "lncRNA" ~ 1,
    TRUE ~ 2
  )
  genes$chrom_rank <- chrom_rank[genes$chrom]
  genes |>
    dplyr::arrange(gene, priority, chrom_rank, start) |>
    dplyr::distinct(gene, .keep_all = TRUE) |>
    dplyr::arrange(chrom_rank, start, end)
}

build_gene_bins <- function(gene_order, bin_size = 80, min_genes = 30) {
  bins <- list()
  for (chrom in unique(gene_order$chrom)) {
    chrom_genes <- gene_order$gene[gene_order$chrom == chrom]
    starts <- seq(1, length(chrom_genes), by = bin_size)
    for (i in seq_along(starts)) {
      genes <- chrom_genes[starts[[i]]:min(starts[[i]] + bin_size - 1, length(chrom_genes))]
      if (length(genes) >= min_genes) {
        bins[[length(bins) + 1]] <- data.frame(
          bin = paste0(chrom, ":", i),
          chrom = chrom,
          bin_no = i,
          gene = genes,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  dplyr::bind_rows(bins)
}

ctc_counts <- get_count_matrix(ctc)
pbmc_counts <- get_count_matrix(pbmc)

gene_order <- parse_gencode_gene_order(PATHS$gene_annotation_gtf)
common_genes <- gene_order |>
  dplyr::filter(gene %in% rownames(ctc_counts), gene %in% rownames(pbmc_counts)) |>
  dplyr::arrange(chrom_rank, start, end)

if (nrow(common_genes) < 1000) {
  stop("Too few common ordered genes for inferCNV-like analysis: ", nrow(common_genes), call. = FALSE)
}

bins_long <- build_gene_bins(common_genes, bin_size = 80, min_genes = 30)
bin_info <- bins_long |>
  dplyr::group_by(bin, chrom, bin_no) |>
  dplyr::summarise(gene_count = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(chrom_rank = match(chrom, c(as.character(1:22), "X", "Y"))) |>
  dplyr::arrange(chrom_rank, bin_no)
write.csv(bin_info, file.path(out_dir, "infercnv_like_gene_bins.csv"), row.names = FALSE)

ordered_genes <- common_genes$gene
gene_index <- match(bins_long$gene, ordered_genes)
bin_index <- match(bins_long$bin, bin_info$bin)
bin_matrix <- Matrix::sparseMatrix(
  i = gene_index,
  j = bin_index,
  x = 1 / bin_info$gene_count[bin_index],
  dims = c(length(ordered_genes), nrow(bin_info)),
  dimnames = list(ordered_genes, bin_info$bin)
)

ctc_norm <- log_normalize_counts(ctc_counts[ordered_genes, , drop = FALSE])
pbmc_norm <- log_normalize_counts(pbmc_counts[ordered_genes, , drop = FALSE])

ctc_bins <- as.matrix(Matrix::t(ctc_norm) %*% bin_matrix)
pbmc_bins <- as.matrix(Matrix::t(pbmc_norm) %*% bin_matrix)

center_by_cell <- function(mat) {
  sweep(mat, 1, rowMeans(mat, na.rm = TRUE), "-")
}

smooth_by_chromosome <- function(mat, bin_info, window = 5) {
  if (window <= 1) return(mat)
  half_window <- floor(window / 2)
  smoothed <- mat
  for (chrom in unique(bin_info$chrom)) {
    idx <- which(bin_info$chrom == chrom)
    if (length(idx) < 2) next
    chrom_mat <- mat[, idx, drop = FALSE]
    chrom_smooth <- chrom_mat
    for (j in seq_along(idx)) {
      local_idx <- max(1, j - half_window):min(length(idx), j + half_window)
      chrom_smooth[, j] <- rowMeans(chrom_mat[, local_idx, drop = FALSE], na.rm = TRUE)
    }
    smoothed[, idx] <- chrom_smooth
  }
  smoothed
}

ref_mean <- colMeans(pbmc_bins, na.rm = TRUE)

# inferCNV-like signal: subtract the normal PBMC reference, remove each cell's
# genome-wide expression offset, smooth locally within chromosomes, then scale
# by a global PBMC-derived background. This avoids imposing similar patterns on
# every chromosome through per-bin standardization.
ctc_signal <- center_by_cell(sweep(ctc_bins, 2, ref_mean, "-"))
pbmc_signal <- center_by_cell(sweep(pbmc_bins, 2, ref_mean, "-"))

ctc_signal <- smooth_by_chromosome(ctc_signal, bin_info, window = 5)
pbmc_signal <- smooth_by_chromosome(pbmc_signal, bin_info, window = 5)

background_scale <- stats::median(abs(as.vector(pbmc_signal)), na.rm = TRUE) / 0.6745
if (!is.finite(background_scale) || background_scale < 0.05) {
  background_scale <- 0.05
}

ctc_z <- ctc_signal / background_scale
pbmc_z <- pbmc_signal / background_scale

ctc_group <- if ("enrichment_method" %in% colnames(ctc@meta.data)) {
  ifelse(tolower(as.character(ctc$enrichment_method)) == "experimental", "Experimental CTC", "Computational CTC")
} else if ("CTC_group" %in% colnames(ctc@meta.data)) {
  as.character(ctc$CTC_group)
} else {
  rep("CTC", ncol(ctc))
}

cnv_scores <- dplyr::bind_rows(
  data.frame(cell = colnames(ctc), group = ctc_group, cnv_proxy_score = rowMeans(abs(ctc_z)), stringsAsFactors = FALSE),
  data.frame(cell = colnames(pbmc), group = "PBMC", cnv_proxy_score = rowMeans(abs(pbmc_z)), stringsAsFactors = FALSE)
)
cnv_scores$group <- factor(cnv_scores$group, levels = c("PBMC", "Experimental CTC", "Computational CTC"))
write.csv(cnv_scores, file.path(out_dir, "cnv_proxy_scores.csv"), row.names = FALSE)

group_matrix <- rbind(
  "PBMC" = colMeans(pbmc_z),
  "Experimental CTC" = colMeans(ctc_z[ctc_group == "Experimental CTC", , drop = FALSE]),
  "Computational CTC" = colMeans(ctc_z[ctc_group == "Computational CTC", , drop = FALSE])
)
write.csv(group_matrix, file.path(out_dir, "infercnv_like_group_mean_matrix.csv"))

set.seed(20260719)
sample_group <- function(mat, n = 120) {
  idx <- sample(seq_len(nrow(mat)), min(n, nrow(mat)))
  if (length(idx) > 2) {
    hc <- hclust(dist(pmin(pmax(mat[idx, , drop = FALSE], -2.5), 2.5)), method = "average")
    idx <- idx[hc$order]
  }
  idx
}

idx_pbmc <- sample_group(pbmc_z, 80)
idx_exp <- sample_group(ctc_z[ctc_group == "Experimental CTC", , drop = FALSE], 130)
idx_comp <- sample_group(ctc_z[ctc_group == "Computational CTC", , drop = FALSE], 130)

exp_global <- which(ctc_group == "Experimental CTC")[idx_exp]
comp_global <- which(ctc_group == "Computational CTC")[idx_comp]

heat_mat <- rbind(
  pbmc_z[idx_pbmc, , drop = FALSE],
  ctc_z[exp_global, , drop = FALSE],
  ctc_z[comp_global, , drop = FALSE]
)
heat_group <- c(rep("PBMC reference", length(idx_pbmc)), rep("Experimental CTC", length(exp_global)), rep("Computational CTC", length(comp_global)))

heat_df <- as.data.frame(as.table(pmin(pmax(heat_mat, -2.5), 2.5)))
colnames(heat_df) <- c("cell_index", "bin", "z")
heat_df$cell_index <- as.integer(heat_df$cell_index)
heat_df$bin_index <- match(heat_df$bin, bin_info$bin)
heat_df$group <- heat_group[heat_df$cell_index]

chrom_bounds <- bin_info |>
  dplyr::mutate(bin_index = dplyr::row_number()) |>
  dplyr::group_by(chrom) |>
  dplyr::summarise(start = min(bin_index), end = max(bin_index), mid = mean(range(bin_index)), .groups = "drop")

p_heat <- ggplot(heat_df, aes(x = bin_index, y = cell_index, fill = z)) +
  geom_tile(width = 1, height = 1) +
  geom_vline(data = chrom_bounds, aes(xintercept = end + 0.5), color = "black", linewidth = 0.12) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
  scale_x_continuous(breaks = chrom_bounds$mid, labels = chrom_bounds$chrom, expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = "Genomic region", y = "Cells", fill = "CNV z", title = "Expression-derived inferCNV-like profile") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave(file.path(out_dir, "infercnv_like_heatmap.pdf"), p_heat, width = 12, height = 7)
ggsave(file.path(out_dir, "infercnv_like_heatmap.jpg"), p_heat, width = 12, height = 7, dpi = 300)

p_cnv <- ggplot(cnv_scores, aes(x = group, y = cnv_proxy_score, fill = group)) +
  geom_violin(width = 0.76, alpha = 0.75, color = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.15, alpha = 0.85) +
  scale_fill_manual(values = c("PBMC" = "#CBD5E1", "Experimental CTC" = "#009E73", "Computational CTC" = "#E69F00")) +
  labs(x = NULL, y = "CNV proxy score", title = "CNV proxy score by evidence group") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(out_dir, "cnv_proxy_score_violin.pdf"), p_cnv, width = 5.2, height = 4.2)
ggsave(file.path(out_dir, "cnv_proxy_score_violin.jpg"), p_cnv, width = 5.2, height = 4.2, dpi = 300)

marker_blocks <- list(
  "Epithelial CTC" = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "KRT7", "ELF3"),
  "EMT/mesenchymal" = c("VIM", "FN1", "ZEB1", "TWIST1", "SNAI1", "SNAI2", "CDH2"),
  "Leukocyte" = c("PTPRC", "CD3D", "CD3E", "LYZ", "MS4A1", "CD68", "NKG7", "LST1")
)

marker_genes <- unique(unlist(marker_blocks))
marker_genes <- marker_genes[marker_genes %in% rownames(ctc_counts) & marker_genes %in% rownames(pbmc_counts)]
gene_to_module <- unlist(lapply(names(marker_blocks), function(module) {
  stats::setNames(rep(module, length(marker_blocks[[module]])), marker_blocks[[module]])
}), use.names = TRUE)
ctc_marker <- as.matrix(log_normalize_counts(ctc_counts[marker_genes, , drop = FALSE]))
pbmc_marker <- as.matrix(log_normalize_counts(pbmc_counts[marker_genes, , drop = FALSE]))

marker_rows <- list()
for (grp in c("Experimental CTC", "Computational CTC")) {
  mat <- ctc_marker[, ctc_group == grp, drop = FALSE]
  for (gene in marker_genes) {
    marker_rows[[length(marker_rows) + 1]] <- data.frame(
      group = grp,
      gene = gene,
      mean_expression = mean(mat[gene, ]),
      pct_expressing = mean(mat[gene, ] > 0) * 100,
      stringsAsFactors = FALSE
    )
  }
}
for (gene in marker_genes) {
  marker_rows[[length(marker_rows) + 1]] <- data.frame(
    group = "PBMC",
    gene = gene,
    mean_expression = mean(pbmc_marker[gene, ]),
    pct_expressing = mean(pbmc_marker[gene, ] > 0) * 100,
    stringsAsFactors = FALSE
  )
}
marker_dot <- dplyr::bind_rows(marker_rows)
marker_dot$module <- unname(gene_to_module[as.character(marker_dot$gene)])
marker_dot$group <- factor(marker_dot$group, levels = c("Experimental CTC", "Computational CTC", "PBMC"))
marker_dot$gene <- factor(marker_dot$gene, levels = marker_genes)
write.csv(marker_dot, file.path(out_dir, "marker_gene_dotplot_data.csv"), row.names = FALSE)

p_dot <- ggplot(marker_dot, aes(x = gene, y = group)) +
  geom_point(aes(size = pct_expressing, color = mean_expression)) +
  facet_grid(. ~ module, scales = "free_x", space = "free_x") +
  scale_size(range = c(1.2, 8.0), name = "% expressing") +
  scale_color_gradientn(colors = c("#DDEAF7", "#74A9CF", "#2B8CBE", "#084081"), name = "Mean expression") +
  labs(x = NULL, y = NULL, title = "Marker-gene expression across confidence groups") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "#F8FAFC", color = "#CBD5E1"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave(file.path(out_dir, "marker_gene_dotplot.pdf"), p_dot, width = 10.5, height = 4.2)
ggsave(file.path(out_dir, "marker_gene_dotplot.jpg"), p_dot, width = 10.5, height = 4.2, dpi = 300)

summary_tbl <- data.frame(
  metric = c("median_cnv_proxy_PBMC", "median_cnv_proxy_Experimental_CTC", "median_cnv_proxy_Computational_CTC"),
  value = c(
    median(cnv_scores$cnv_proxy_score[cnv_scores$group == "PBMC"], na.rm = TRUE),
    median(cnv_scores$cnv_proxy_score[cnv_scores$group == "Experimental CTC"], na.rm = TRUE),
    median(cnv_scores$cnv_proxy_score[cnv_scores$group == "Computational CTC"], na.rm = TRUE)
  )
)
write.csv(summary_tbl, file.path(out_dir, "confidence_validation_summary.csv"), row.names = FALSE)

# Optional official infercnv package workflow.
# Set RUN_OFFICIAL_INFERCNV <- TRUE only when infercnv is installed and enough
# memory is available. The lightweight inferCNV-like output above is used for
# database confidence figures by default.
RUN_OFFICIAL_INFERCNV <- FALSE
if (RUN_OFFICIAL_INFERCNV) {
  required_packages(c("infercnv"))
  infercnv_dir <- file.path(out_dir, "official_infercnv")
  dir.create(infercnv_dir, showWarnings = FALSE, recursive = TRUE)

  common_for_infercnv <- intersect(rownames(ctc_counts), rownames(pbmc_counts))
  count_all <- Matrix::cbind2(ctc_counts[common_for_infercnv, ], pbmc_counts[common_for_infercnv, ])
  annotations <- data.frame(
    cell = colnames(count_all),
    group = c(ctc_group, rep("PBMC", ncol(pbmc_counts))),
    stringsAsFactors = FALSE
  )
  annotation_file <- file.path(infercnv_dir, "cell_annotations.tsv")
  write.table(annotations, annotation_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  gene_order_file <- file.path(infercnv_dir, "gene_order.tsv")
  gene_order_tbl <- common_genes |>
    dplyr::filter(gene %in% common_for_infercnv) |>
    dplyr::select(gene, chrom, start, end)
  write.table(gene_order_tbl, gene_order_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  infer_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = count_all,
    annotations_file = annotation_file,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = c("PBMC")
  )
  infercnv::run(
    infer_obj,
    cutoff = 0.1,
    out_dir = infercnv_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    num_threads = max(1, parallel::detectCores() - 1)
  )
}

message("Confidence validation outputs saved to: ", out_dir)
