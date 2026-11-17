#!/usr/bin/env Rscript

# Single-cell cancer-vs-all-other-cancers DEG and GO/KEGG analysis.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(normalizePath(file.path(script_dir, "..", "01_config", "01_setup_paths.R"), mustWork = TRUE))
required_packages(c("Seurat", "ggplot2", "dplyr", "ggrepel", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))

scored_rds <- file.path(PATHS$output, "03_SEED_scRNA_scored.rds")
if (!file.exists(scored_rds)) {
  scored_rds <- file.path(PATHS$output, "04_SEED_scRNA_scored.rds")
}
obj <- require_rds(scored_rds, "scored SEED scRNA object")
obj <- set_default_assay(obj)
meta <- obj@meta.data
out_dir <- PATHS$sc_cancer_output

symbol2entrez <- function(genes) {
  if (length(genes) < 1) return(character(0))
  suppressMessages({
    mp <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  })
  unique(mp$ENTREZID)
}

run_enrich <- function(genes, prefix, out_dir) {
  eg <- symbol2entrez(genes)
  if (length(eg) < 5) return(invisible(NULL))
  go <- enrichGO(gene = eg, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  kegg <- enrichKEGG(gene = eg, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  if (!is.null(go) && nrow(go) > 0) {
    write.csv(as.data.frame(go), file.path(out_dir, paste0(prefix, "_GO_BP.csv")), row.names = FALSE)
    if (nrow(go) > 1) save_plot(dotplot(go, showCategory = 15, title = paste0(prefix, " GO BP")), file.path(out_dir, paste0(prefix, "_GO_dotplot")), 10, 8)
  }
  if (!is.null(kegg) && nrow(kegg) > 0) {
    write.csv(as.data.frame(kegg), file.path(out_dir, paste0(prefix, "_KEGG.csv")), row.names = FALSE)
    if (nrow(kegg) > 1) save_plot(dotplot(kegg, showCategory = 15, title = paste0(prefix, " KEGG")), file.path(out_dir, paste0(prefix, "_KEGG_dotplot")), 10, 8)
  }
}

plot_volcano <- function(deg, prefix, title) {
  dat <- deg |>
    dplyr::filter(!is.na(avg_log2FC), !is.na(p_val_adj)) |>
    dplyr::mutate(
      negLogFDR = -log10(pmax(p_val_adj, .Machine$double.xmin)),
      group = dplyr::case_when(
        avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Up",
        avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Down",
        TRUE ~ "NS"
      )
    )
  lab <- dplyr::bind_rows(
    dat |> dplyr::filter(group == "Up") |> dplyr::arrange(p_val_adj) |> head(10),
    dat |> dplyr::filter(group == "Down") |> dplyr::arrange(p_val_adj) |> head(10)
  )
  p <- ggplot(dat, aes(avg_log2FC, negLogFDR)) +
    geom_point(aes(color = group), alpha = 0.7, size = 1.4) +
    scale_color_manual(values = c(Up = "#D73027", Down = "#4575B4", NS = "grey80")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3.4, max.overlaps = Inf) +
    labs(x = "log2 fold-change", y = "-log10 adjusted P", title = title, color = NULL) +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p, file.path(out_dir, prefix), 9, 7)
}

min_cells_per_group <- 50
max_size_ratio <- 2

for (target in na.omit(unique(as.character(meta$cancer)))) {
  message("Cancer-vs-others: ", target)
  cells_a <- rownames(meta)[meta$cancer == target]
  cells_o <- rownames(meta)[meta$cancer != target & !is.na(meta$cancer)]
  if (length(cells_a) < min_cells_per_group || length(cells_o) < min_cells_per_group) next
  if (max(length(cells_a), length(cells_o)) / min(length(cells_a), length(cells_o)) > max_size_ratio) {
    set.seed(123)
    n <- min(length(cells_a), length(cells_o))
    cells_a <- sample(cells_a, n)
    cells_o <- sample(cells_o, n)
  }
  cells_use <- c(cells_a, cells_o)
  sub_obj <- subset(obj, cells = cells_use)
  sub_obj$temp_cancer_group <- ifelse(colnames(sub_obj) %in% cells_a, target, "Others")

  deg <- FindMarkers(sub_obj, ident.1 = target, ident.2 = "Others", group.by = "temp_cancer_group", test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, only.pos = FALSE)
  if (nrow(deg) == 0) next

  deg$gene <- rownames(deg)
  deg <- deg[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
  sig <- subset(deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
  prefix <- paste0(target, "_vs_AllOtherCancers")
  write.csv(deg, file.path(out_dir, paste0(prefix, "_DEGs_all.csv")), row.names = FALSE)
  write.csv(sig, file.path(out_dir, paste0(prefix, "_DEGs_significant.csv")), row.names = FALSE)
  plot_volcano(deg, paste0(prefix, "_volcano"), paste0(target, " vs all other cancers"))
  if (nrow(sig) >= 5) {
    run_enrich(sig$gene[sig$avg_log2FC > 0], paste0(prefix, "_up"), out_dir)
    run_enrich(sig$gene[sig$avg_log2FC < 0], paste0(prefix, "_down"), out_dir)
    run_enrich(sig$gene, paste0(prefix, "_all"), out_dir)
  }
}

message("Cancer-vs-others scRNA analysis finished: ", out_dir)
