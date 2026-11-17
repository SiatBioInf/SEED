#!/usr/bin/env Rscript

# State-level High/Low DEG and GO/KEGG analysis for SEED scRNA CTC data.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(normalizePath(file.path(script_dir, "..", "01_config", "01_setup_paths.R"), mustWork = TRUE))
required_packages(c("Seurat", "ggplot2", "dplyr", "tidyr", "ggrepel", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))

scored_rds <- file.path(PATHS$output, "03_SEED_scRNA_scored.rds")
if (!file.exists(scored_rds)) {
  scored_rds <- file.path(PATHS$output, "04_SEED_scRNA_scored.rds")
}
obj <- require_rds(scored_rds, "scored SEED scRNA object")
obj <- set_default_assay(obj)
meta <- obj@meta.data
score_cols <- intersect(STATE_THEMES$score, grep("_Score$", colnames(meta), value = TRUE))
out_dir <- PATHS$sc_state_output

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

for (score_name in score_cols) {
  message("Analyzing ", score_name)
  group_col <- paste0(score_name, "_Group")
  if (!group_col %in% colnames(meta)) next

  write.csv(table(meta[[group_col]]), file.path(out_dir, paste0(score_name, "_cell_count_simple.csv")))
  if ("cancer" %in% colnames(meta)) {
    write.csv(table(meta[[group_col]], meta$cancer), file.path(out_dir, paste0(score_name, "_cell_count_by_cancer.csv")))
  }

  p_vln <- VlnPlot(obj, features = score_name, group.by = group_col, pt.size = 0) +
    ggtitle(paste0(score_name, " High/Low"))
  save_plot(p_vln, file.path(out_dir, paste0(score_name, "_violin")), 8, 6)

  if ("umap" %in% names(obj@reductions)) {
    p_umap <- DimPlot(obj, reduction = "umap", group.by = group_col) + ggtitle(paste0(score_name, " High/Low UMAP"))
    save_plot(p_umap, file.path(out_dir, paste0(score_name, "_umap")), 9, 6)
  }

  deg <- FindMarkers(obj, ident.1 = "High", ident.2 = "Low", group.by = group_col, test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, only.pos = FALSE)
  if (nrow(deg) == 0) next
  deg$gene <- rownames(deg)
  deg <- deg[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
  sig <- subset(deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
  write.csv(deg, file.path(out_dir, paste0(score_name, "_DEGs_all.csv")), row.names = FALSE)
  write.csv(sig, file.path(out_dir, paste0(score_name, "_DEGs_significant.csv")), row.names = FALSE)
  plot_volcano(deg, paste0(score_name, "_volcano"), paste0(score_name, " High vs Low"))

  if (nrow(sig) >= 5) {
    run_enrich(sig$gene[sig$avg_log2FC > 0], paste0(score_name, "_up"), out_dir)
    run_enrich(sig$gene[sig$avg_log2FC < 0], paste0(score_name, "_down"), out_dir)
    run_enrich(sig$gene, paste0(score_name, "_all"), out_dir)
  }

  if ("cancer" %in% colnames(meta)) {
    for (cancer in na.omit(unique(as.character(meta$cancer)))) {
      cells_cancer <- rownames(meta)[meta$cancer == cancer]
      high_cells <- intersect(cells_cancer, rownames(meta)[meta[[group_col]] == "High"])
      low_cells <- intersect(cells_cancer, rownames(meta)[meta[[group_col]] == "Low"])
      if (length(high_cells) < min_cells_per_group || length(low_cells) < min_cells_per_group) next
      if (max(length(high_cells), length(low_cells)) / min(length(high_cells), length(low_cells)) > max_size_ratio) {
        set.seed(123)
        n <- min(length(high_cells), length(low_cells))
        high_cells <- sample(high_cells, n)
        low_cells <- sample(low_cells, n)
      }
      cells_use <- c(high_cells, low_cells)
      sub_obj <- subset(obj, cells = cells_use)
      deg_c <- FindMarkers(sub_obj, ident.1 = "High", ident.2 = "Low", group.by = group_col, test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, only.pos = FALSE)
      if (nrow(deg_c) == 0) next
      deg_c$gene <- rownames(deg_c)
      deg_c <- deg_c[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      sig_c <- subset(deg_c, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
      prefix <- paste0(score_name, "_", cancer)
      write.csv(deg_c, file.path(out_dir, paste0(prefix, "_DEGs_all.csv")), row.names = FALSE)
      write.csv(sig_c, file.path(out_dir, paste0(prefix, "_DEGs_significant.csv")), row.names = FALSE)
      plot_volcano(deg_c, paste0(prefix, "_volcano"), paste0(score_name, " High vs Low in ", cancer))
      if (nrow(sig_c) >= 5) {
        run_enrich(sig_c$gene[sig_c$avg_log2FC > 0], paste0(prefix, "_up"), out_dir)
        run_enrich(sig_c$gene[sig_c$avg_log2FC < 0], paste0(prefix, "_down"), out_dir)
        run_enrich(sig_c$gene, paste0(prefix, "_all"), out_dir)
      }
    }
  }
}

saveRDS(obj, file.path(PATHS$output, "05_SEED_scRNA_state_groups_and_results.rds"))
message("State analysis finished: ", out_dir)
