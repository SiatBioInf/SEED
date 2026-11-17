#!/usr/bin/env Rscript

# Bulk RNA cancer-vs-all-other-cancers DEG and GO/KEGG analysis.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(normalizePath(file.path(script_dir, "..", "01_config", "01_setup_paths.R"), mustWork = TRUE))
required_packages(c("Seurat", "DESeq2", "ggplot2", "dplyr", "ggrepel", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))

bulk <- require_rds(PATHS$bulk_rds, "SEED bulk RNA")
bulk <- set_default_assay(bulk)
meta <- as.data.frame(bulk@meta.data)
bulk_assay <- DefaultAssay(bulk)
counts <- tryCatch(
  GetAssayData(bulk, assay = bulk_assay, layer = "counts"),
  error = function(e) GetAssayData(bulk, assay = bulk_assay, slot = "counts")
)
counts <- round(counts[, rownames(meta), drop = FALSE])
out_dir <- PATHS$bulk_cancer_output

symbol2entrez <- function(genes) {
  if (length(genes) < 1) return(character(0))
  suppressMessages({
    mp <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  })
  unique(mp$ENTREZID)
}

run_bulk_enrich <- function(genes, prefix) {
  eg <- symbol2entrez(genes)
  if (length(eg) < 5) return(invisible(NULL))
  go <- enrichGO(gene = eg, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  kegg <- enrichKEGG(gene = eg, organism = "hsa", pvalueCutoff = 0.05)
  if (!is.null(go) && nrow(go) > 0) {
    write.csv(as.data.frame(go), file.path(out_dir, paste0(prefix, "_GO_BP.csv")), row.names = FALSE)
    save_plot(dotplot(go, showCategory = 15, title = paste0(prefix, " GO BP")), file.path(out_dir, paste0(prefix, "_GO_dotplot")), 10, 8)
  }
  if (!is.null(kegg) && nrow(kegg) > 0) {
    write.csv(as.data.frame(kegg), file.path(out_dir, paste0(prefix, "_KEGG.csv")), row.names = FALSE)
    save_plot(dotplot(kegg, showCategory = 15, title = paste0(prefix, " KEGG")), file.path(out_dir, paste0(prefix, "_KEGG_dotplot")), 10, 8)
  }
}

plot_bulk_volcano <- function(de, prefix, title) {
  dat <- de |>
    dplyr::filter(!is.na(log2FoldChange), !is.na(padj)) |>
    dplyr::mutate(
      negLogFDR = -log10(pmax(padj, .Machine$double.xmin)),
      group = dplyr::case_when(
        log2FoldChange > 0.5 & padj < 0.05 ~ "Up",
        log2FoldChange < -0.5 & padj < 0.05 ~ "Down",
        TRUE ~ "NS"
      )
    )
  lab <- dplyr::bind_rows(
    dat |> dplyr::filter(group == "Up") |> dplyr::arrange(padj) |> head(10),
    dat |> dplyr::filter(group == "Down") |> dplyr::arrange(padj) |> head(10)
  )
  p <- ggplot(dat, aes(log2FoldChange, negLogFDR)) +
    geom_point(aes(color = group), alpha = 0.7, size = 1.3) +
    scale_color_manual(values = c(Up = "#D73027", Down = "#4575B4", NS = "grey80")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3.4, max.overlaps = Inf) +
    labs(x = "log2 fold-change", y = "-log10 FDR", title = title, color = NULL) +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p, file.path(out_dir, prefix), 9, 7)
}

min_sample_per_group <- 3
max_ratio <- 3

write.csv(dplyr::count(meta, cancer, name = "sample_n"), file.path(out_dir, "bulkRNA_cancer_sample_count.csv"), row.names = FALSE)

for (target in na.omit(unique(as.character(meta$cancer)))) {
  message("Bulk cancer-vs-others: ", target)
  samp_a <- rownames(meta)[meta$cancer == target]
  samp_o <- rownames(meta)[meta$cancer != target & !is.na(meta$cancer)]
  if (length(samp_a) < min_sample_per_group || length(samp_o) < min_sample_per_group) next
  if (max(length(samp_a), length(samp_o)) / min(length(samp_a), length(samp_o)) > max_ratio) {
    set.seed(123)
    n <- min(length(samp_a), length(samp_o))
    samp_a <- sample(samp_a, n)
    samp_o <- sample(samp_o, n)
  }

  samp_use <- c(samp_a, samp_o)
  sub_meta <- meta[samp_use, , drop = FALSE]
  sub_meta$compare_group <- factor(ifelse(sub_meta$cancer == target, target, "Other"), levels = c("Other", target))
  dds <- DESeqDataSetFromMatrix(countData = counts[, samp_use, drop = FALSE], colData = sub_meta, design = ~ compare_group)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("compare_group", target, "Other"), lfcThreshold = 0.25)
  de <- as.data.frame(res) |>
    dplyr::mutate(gene = rownames(.)) |>
    dplyr::select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
  sig <- dplyr::filter(de, padj < 0.05 & abs(log2FoldChange) > 0.5)

  prefix <- paste0(target, "_vs_AllOthers")
  write.csv(de, file.path(out_dir, paste0(prefix, "_DE_all.csv")), row.names = FALSE)
  write.csv(sig, file.path(out_dir, paste0(prefix, "_DE_sig.csv")), row.names = FALSE)
  plot_bulk_volcano(de, paste0(prefix, "_volcano"), paste0(target, " vs all other cancers"))
  run_bulk_enrich(sig$gene[sig$log2FoldChange > 0], paste0(prefix, "_up"))
  run_bulk_enrich(sig$gene[sig$log2FoldChange < 0], paste0(prefix, "_down"))
  run_bulk_enrich(sig$gene, paste0(prefix, "_all_sig"))
}

message("Bulk analysis finished: ", out_dir)
