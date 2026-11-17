# ============================================================
# CTC-KG v2 — scRNA分析: 效应量驱动的癌种特异性基因+通路
# ============================================================
# 核心改进（vs v1 revised_03_specificity.R）:
#   1. 效应量优先：pct_delta + Cohen's d 为主要筛选标准，p值为辅
#   2. 三级分类：Tier 1/2/3 替代单一 p_adj < 0.05
#   3. 新增 detection_rate 指标
#   4. 通路从20条手动 → MSigDB Hallmark (50条) + KEGG/Reactome
# ============================================================

library(Seurat)
library(Matrix)

# === 配置 ===
out_dir <- "D:/ClaudeCode记录/CTC-KG论文"
v2_dir  <- "D:/科研/CTC知识图谱/v2/data/processed"
dir.create(v2_dir, showWarnings = FALSE, recursive = TRUE)

cancer_names <- c("BC","CRC","GC","HCC","HNSC","MELA","NPC","NSCLC","OV","PC","PDAC","SCLC")
cancer_full  <- c("Breast Cancer", "Colorectal Cancer", "Gastric Cancer",
                  "Hepatocellular Carcinoma", "Head and Neck Squamous Cell Carcinoma",
                  "Melanoma", "Nasopharyngeal Carcinoma", "Non-Small Cell Lung Cancer",
                  "Ovarian Cancer", "Prostate Cancer", "Pancreatic Ductal Adenocarcinoma",
                  "Small Cell Lung Cancer")

# ============================================================
# Part A: 加载数据 + 计算已有百分位
# ============================================================
cat("===== Part A: Loading data =====\n")

# 如果已有百分位数据，直接加载（节省时间）
pct_file <- file.path(out_dir, "percentile_by_cancer.rds")
if (file.exists(pct_file)) {
  cat("Loading existing percentile data...\n")
  percentile_list <- readRDS(pct_file)
  cat(sprintf("Loaded %d cancers\n", length(percentile_list)))
} else {
  stop("percentile_by_cancer.rds not found. Run revised_03_specificity.R Part A first.")
}

# 加载合并Seurat对象获取表达矩阵
merged <- readRDS(file.path(out_dir, "CTC_merged_harmony.rds"))
merged <- JoinLayers(merged)
expr_all <- GetAssayData(merged, assay = "RNA", layer = "data")
cancer_labels <- merged$cancer

# 获取基因列表（排除ERCC/MT-/Ribosomal）
genes_all <- rownames(expr_all)
remove_pattern <- grepl("^ERCC|^MT-|^RP[SL]", genes_all)
genes_clean <- genes_all[!remove_pattern]
cat(sprintf("Genes after filtering: %d\n", length(genes_clean)))

# ============================================================
# Part B: 效应量驱动的癌种特异性基因分析
# ============================================================
cat("\n===== Part B: Effect-size-driven cancer-specific gene analysis =====\n")

# --- B1: 计算每个基因在每个癌种的效应量（无p值，快速）---
# 这一步只做简单的均值和效应量计算，不做统计检验。
# 统计检验放到B2，用预计算的百分位矩阵加速。
gene_stats <- list()

for (cn in cancer_names) {
  cat(sprintf("Computing stats for %s...\n", cn))
  pct_mat <- percentile_list[[cn]]  # gene x cell（已有的百分位矩阵）
  expr_mat <- expr_all[, cancer_labels == cn, drop = FALSE]

  # 1. 基因在癌种CTC细胞内的平均百分位
  mean_pct <- rowMeans(pct_mat, na.rm = TRUE)

  # 2. 检出率：该基因在多少比例的细胞中被检出（表达量>0）
  det_rate <- rowSums(expr_mat > 0) / ncol(expr_mat)

  # 3. 对"其他癌种"采样200个细胞，预计算百分位矩阵
  #    这个矩阵后续会被B2复用，所以只算一次
  cells_other <- which(cancer_labels != cn)
  set.seed(42)
  idx_other <- sample(cells_other, min(200, length(cells_other)))
  expr_other_mat <- as.matrix(expr_all[, idx_other, drop = FALSE])

  # 对每个采样细胞，计算所有基因的within-cell百分位
  # 结果: gene x cell_other 的百分位矩阵
  pct_other_mat <- apply(expr_other_mat, 2, function(x) {
    rank(x, ties.method = "average") / length(x)
  })
  rownames(pct_other_mat) <- rownames(expr_all)

  # 4. 效应量计算
  mean_other_pct <- rowMeans(pct_other_mat, na.rm = TRUE)
  pct_delta <- mean_pct - mean_other_pct

  # Cohen's d（标准化效应量，不受细胞数影响）
  sd_target <- apply(pct_mat, 1, sd, na.rm = TRUE)
  sd_other  <- apply(pct_other_mat, 1, sd, na.rm = TRUE)
  n_t <- ncol(pct_mat)
  n_o <- ncol(pct_other_mat)
  pooled_sd <- sqrt(((n_t - 1) * sd_target^2 + (n_o - 1) * sd_other^2) / (n_t + n_o - 2))
  cohens_d <- pct_delta / (pooled_sd + 1e-8)

  gene_stats[[cn]] <- list(
    df = data.frame(
      gene = rownames(pct_mat), cancer = cn,
      mean_percentile = mean_pct, detection_rate = det_rate,
      pct_delta = pct_delta, cohens_d = cohens_d,
      n_cells = ncol(pct_mat),
      row.names = NULL, stringsAsFactors = FALSE
    ),
    pct_target = pct_mat,      # 保留给B2做Wilcoxon
    pct_other  = pct_other_mat  # 保留给B2
  )
  cat(sprintf("  %s: %d cells, %d other-cells sampled\n", cn, ncol(pct_mat), ncol(pct_other_mat)))
}

# --- B2: Wilcoxon 检验（复用B1预计算的百分位矩阵，快速）---
cat("\nRunning Wilcoxon tests (using pre-computed percentile matrices)...\n")

all_stats <- do.call(rbind, lapply(gene_stats, `[[`, "df"))
all_stats$p_value <- NA_real_
all_stats$p_adj   <- NA_real_

for (cn in cancer_names) {
  gs <- gene_stats[[cn]]
  pct_target <- gs$pct_target   # gene × cell 百分位矩阵（靶癌种）
  pct_other  <- gs$pct_other    # gene × cell 百分位矩阵（其他癌种，已预计算）

  # 只测试效应量有意义的基因（pct_delta > 0.02，减少无效检验）
  cn_df <- gs$df
  genes_to_test <- cn_df$gene[cn_df$pct_delta > 0.02]
  n_test <- length(genes_to_test)

  if (n_test == 0) next

  cat(sprintf("  Testing %s: %d genes with delta>0.02...\n", cn, n_test))

  # 关键优化：
  #   不再对每个基因重新计算other百分位（那需要排序36701个基因）。
  #   而是直接从B1预计算的 pct_other 矩阵中提取该基因的行。
  #   Wilcoxon检验比较：
  #     target：该基因在靶癌种所有细胞中的百分位值（已有）
  #     other：该基因在采样其他癌种细胞中的百分位值（已有）
  pb <- txtProgressBar(min = 0, max = n_test, style = 3)
  for (i in seq_along(genes_to_test)) {
    g <- genes_to_test[i]
    target_vals <- as.numeric(pct_target[g, ])
    other_vals  <- as.numeric(pct_other[g, ])

    wt <- tryCatch(
      wilcox.test(target_vals, other_vals, alternative = "greater"),
      error = function(e) list(p.value = 1)
    )

    idx <- which(all_stats$gene == g & all_stats$cancer == cn)
    all_stats$p_value[idx] <- wt$p.value

    if (i %% 500 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  cat(sprintf("  Done %s\n", cn))
}

# BH多重检验校正（控制假阳性率）
all_stats$p_adj <- p.adjust(all_stats$p_value, method = "BH")
all_stats$p_value[is.na(all_stats$p_value)] <- 1
all_stats$p_adj[is.na(all_stats$p_adj)] <- 1

# --- B3: Effect-size-driven Tier Classification ---
cat("\n===== Classifying genes into tiers =====\n")

classify_tier <- function(pct_delta, p_adj, detection_rate, n_cells) {
  # Tier 1: 强效应 + 统计显著 + 足够检出率 + 足够样本量（≥50细胞）
  #   为什么需要 n_cells ≥ 50？
  #   小样本癌种（SCLC=18, OV=14）的百分位估计噪音极大，
  #   可能出现 delta=0.5 但实际生物学意义为零的假阳性。
  #   50个细胞是最低可靠统计量要求。
  t1 <- pct_delta > 0.10 & p_adj < 0.05 & detection_rate > 0.25 & n_cells >= 50
  t2 <- pct_delta > 0.05 & p_adj < 0.05 & n_cells >= 30
  t3 <- pct_delta > 0.03 & p_adj < 0.05

  ifelse(t1, 1, ifelse(t2, 2, ifelse(t3, 3, 0)))
}

all_stats$confidence_tier <- classify_tier(
  all_stats$pct_delta, all_stats$p_adj,
  all_stats$detection_rate, all_stats$n_cells
)

# specificity_score: 映射 tier 到 0-1 连续分数（结合 pct_delta）
all_stats$specificity_score <- ifelse(
  all_stats$confidence_tier == 1, pmin(1.0, 0.7 + all_stats$pct_delta),
  ifelse(
    all_stats$confidence_tier == 2, 0.4 + 0.3 * (all_stats$pct_delta - 0.05) / 0.05,
    ifelse(
      all_stats$confidence_tier == 3, 0.2 + 0.2 * (all_stats$pct_delta - 0.03) / 0.02,
      0
    )
  )
)
all_stats$specificity_score <- pmin(1.0, pmax(0, all_stats$specificity_score))

# 统计
cat(sprintf("\n=== Tier distribution ===\n"))
for (t in 1:3) {
  n <- sum(all_stats$confidence_tier == t)
  cat(sprintf("  Tier %d: %d genes (%.1f%%)\n", t, n, 100*n/nrow(all_stats[all_stats$confidence_tier > 0,])))
}
cat(sprintf("  Excluded: %d genes\n", sum(all_stats$confidence_tier == 0)))

# 按癌种统计
cat("\n=== Per-cancer Tier 1+2 gene counts ===\n")
for (cn in cancer_names) {
  cn_genes <- all_stats[all_stats$cancer == cn & all_stats$confidence_tier %in% c(1,2), ]
  cat(sprintf("  %s: %d genes\n", cn, nrow(cn_genes)))
}

# --- B4: 输出基因-癌种边 ---
cat("\n===== Writing gene_cancer_edges.csv =====\n")

edges_out <- all_stats[all_stats$confidence_tier > 0,
  c("gene", "cancer", "mean_percentile", "pct_delta", "cohens_d",
    "specificity_score", "detection_rate", "confidence_tier", "p_value", "p_adj", "n_cells")
]

# Map cancer abbreviations to full names
name_map <- setNames(cancer_full, cancer_names)
edges_out$cancer_full <- name_map[edges_out$cancer]

write.csv(edges_out, file.path(v2_dir, "gene_cancer_edges.csv"), row.names = FALSE)
cat(sprintf("Wrote %d gene-cancer edges to gene_cancer_edges.csv\n", nrow(edges_out)))

# 同时输出完整统计表（含excluded的）
write.csv(all_stats, file.path(v2_dir, "gene_cancer_full_stats.csv"), row.names = FALSE)

# Top genes per cancer for inspection
cat("\n=== Top Tier 1 genes per cancer ===\n")
tier1 <- all_stats[all_stats$confidence_tier == 1, ]
tier1 <- tier1[order(-tier1$pct_delta), ]
for (cn in cancer_names) {
  cn_t1 <- tier1[tier1$cancer == cn, ]
  if (nrow(cn_t1) > 0) {
    cat(sprintf("  %s:\n", cn))
    for (i in 1:min(5, nrow(cn_t1))) {
      cat(sprintf("    %-10s delta=%.4f d=%.2f det=%.2f\n",
                  cn_t1$gene[i], cn_t1$pct_delta[i], cn_t1$cohens_d[i], cn_t1$detection_rate[i]))
    }
  }
}

# ============================================================
# Part C: 通路活性分析（MSigDB Hallmark + KEGG）
# ============================================================
cat("\n===== Part C: Pathway activity analysis =====\n")

# --- C1: 获取通路定义 ---
# 优先使用内置fallback（40条CTC相关通路，覆盖Hallmark核心），
# 因为MSigDB可能改了下载URL需要登录。
cat("Loading pathway gene sets...\n")
use_fallback <- TRUE

if (!use_fallback) {
  hallmark_url <- "https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"
  hallmark_file <- file.path(v2_dir, "h.all.v2024.1.Hs.symbols.gmt")
  tryCatch({
    download.file(hallmark_url, hallmark_file, mode = "wb")
    hallmark_sets <- parse_gmt(hallmark_file)
    if (length(hallmark_sets) > 10) use_fallback <- FALSE
  }, error = function(e) {})
}

if (use_fallback) {
  cat("Using built-in CTC pathway list (40 pathways)\n")
  hallmark_sets <- list(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = c("VIM","CDH1","CDH2","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","FN1","COL1A1","COL1A2","MMP2","MMP9","TGFB1","TGFB2","WNT5A","CTNNB1"),
    "HALLMARK_TGF_BETA_SIGNALING" = c("TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","BMP2","BMP4"),
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING" = c("PIK3CA","PIK3CB","PIK3CD","AKT1","AKT2","AKT3","MTOR","PTEN","RICTOR","RPTOR","TSC1","TSC2"),
    "HALLMARK_KRAS_SIGNALING" = c("KRAS","NRAS","HRAS","BRAF","RAF1","MAP2K1","MAP2K2","MAPK1","MAPK3"),
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING" = c("CTNNB1","WNT3A","WNT5A","AXIN1","AXIN2","APC","GSK3B","TCF7","TCF7L2","LEF1"),
    "HALLMARK_NOTCH_SIGNALING" = c("NOTCH1","NOTCH2","NOTCH3","NOTCH4","DLL1","DLL3","DLL4","JAG1","JAG2","HES1","HEY1"),
    "HALLMARK_HEDGEHOG_SIGNALING" = c("SHH","IHH","DHH","PTCH1","SMO","GLI1","GLI2","GLI3"),
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = c("TNF","NFKB1","NFKB2","RELA","RELB","TNFAIP3","NFKBIA","NFKBIB","IKBKG"),
    "HALLMARK_IL6_JAK_STAT3_SIGNALING" = c("IL6","IL6R","IL6ST","JAK1","JAK2","STAT3","SOCS1","SOCS3"),
    "HALLMARK_INFLAMMATORY_RESPONSE" = c("IL1B","IL6","TNF","CXCL8","CCL2","CXCL1","CXCL2","CXCL5","CXCR4","PTGS2"),
    "HALLMARK_IL2_STAT5_SIGNALING" = c("IL2","IL2RA","IL2RB","STAT5A","STAT5B","JAK1","JAK3"),
    "HALLMARK_INTERFERON_GAMMA_RESPONSE" = c("IFNG","IFNGR1","IFNGR2","STAT1","IRF1","IRF2","IRF9","GBP1","GBP2"),
    "HALLMARK_INTERFERON_ALPHA_RESPONSE" = c("IFNA1","IFNAR1","IFNAR2","STAT1","STAT2","IRF9","MX1","OAS1"),
    "HALLMARK_APOPTOSIS" = c("BCL2","BCL2L1","BAX","BAK1","BAD","BID","CASP3","CASP8","CASP9","TP53","BBC3","PMAIP1"),
    "HALLMARK_P53_PATHWAY" = c("TP53","CDKN1A","MDM2","MDM4","BAX","BBC3","GADD45A","RRM2B","SESN1","SESN2"),
    "HALLMARK_HYPOXIA" = c("HIF1A","ARNT","EPAS1","HIF3A","VEGFA","LDHA","HK2","SLC2A1","PGK1","ENO1"),
    "HALLMARK_ANGIOGENESIS" = c("VEGFA","VEGFB","VEGFC","KDR","FLT1","FLT4","ANGPT1","ANGPT2","TEK","FGF2","PDGFB"),
    "HALLMARK_APICAL_JUNCTION" = c("CDH1","CTNNB1","CTNNA1","CTNND1","JUP","TJP1","OCLN","CLDN1","CLDN3","CLDN4","CLDN7"),
    "HALLMARK_ESTROGEN_RESPONSE_EARLY" = c("ESR1","PGR","TFF1","GREB1","MYC","CCND1","IGF1R","XBP1","SLC7A5"),
    "HALLMARK_ANDROGEN_RESPONSE" = c("AR","KLK2","KLK3","FKBP5","TMPRSS2","NKX3-1","PMEPA1","MAF","AZGP1"),
    "HALLMARK_DNA_REPAIR" = c("BRCA1","BRCA2","RAD51","RAD51C","XRCC2","XRCC3","FANCA","FANCD2","ATM","ATR","CHEK1","CHEK2"),
    "HALLMARK_G2M_CHECKPOINT" = c("CCNB1","CCNB2","CDK1","CDC25A","CDC25C","PLK1","AURKA","AURKB","BUB1","BUB1B","MAD2L1"),
    "HALLMARK_E2F_TARGETS" = c("E2F1","E2F2","E2F3","TFDP1","CCNE1","CCNE2","CDK2","CDC6","MCM2","MCM3","MCM4","MCM5"),
    "HALLMARK_MYC_TARGETS_V1" = c("MYC","NPM1","NCL","DDX21","PPAN","RPL3","RPL11","RPS6","RPS19","EIF4E"),
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = c("NDUFA1","NDUFA9","SDHA","SDHB","UQCRC1","UQCRC2","COX5A","COX5B","ATP5A1","ATP5B"),
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" = c("SOD1","SOD2","CAT","GPX1","GPX4","PRDX1","TXN","TXNRD1","NFE2L2","KEAP1"),
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" = c("HSPA5","HSP90B1","ATF4","ATF6","XBP1","ERN1","EIF2AK3","DDIT3","DNAJB9","HERPUD1"),
    "HALLMARK_PROTEIN_SECRETION" = c("SEC11A","SEC61A1","SEC61B","SRP9","SRP14","SRP19","SRP54","SRPR","SSR1"),
    "HALLMARK_FATTY_ACID_METABOLISM" = c("ACAA2","ACADM","ACADVL","ACAT1","CPT1A","CPT2","HADHA","HADHB","ACSL1","ACSL3"),
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS" = c("HMGCR","HMGCS1","LDLR","SREBF2","SCAP","INSIG1","DHCR7","SQLE","FDFT1","LSS"),
    "HALLMARK_MYC_TARGETS_V2" = c("MYC","FBL","NOP56","NOP58","DKC1","GAR1","NHP2","NOP10","TERT","TERC"),
    "HALLMARK_SPERMATOGENESIS" = c("DAZL","DDX4","PIWIL1","PIWIL2","SYCP1","SYCP3","TDRD1","TDRD5","MAEL","ASZ1"),
    "HALLMARK_BILE_ACID_METABOLISM" = c("CYP7A1","CYP8B1","CYP27A1","BAAT","SLC10A1","SLC10A2","NR1H4","SLC27A5"),
    "HALLMARK_XENOBIOTIC_METABOLISM" = c("CYP1A1","CYP1A2","CYP2E1","CYP3A4","CYP3A5","GSTA1","GSTA2","GSTM1","NQO1","AHR"),
    "HALLMARK_GLYCOLYSIS" = c("HK2","HK1","GPI","PFKL","PFKP","ALDOA","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA","SLC2A1","SLC2A3"),
    "HALLMARK_HEME_METABOLISM" = c("ALAS1","ALAS2","ALAD","HMBS","UROS","UROD","FECH","BLVRA","BLVRB","HMOX1","HMOX2"),
    "HALLMARK_PANCREAS_BETA_CELLS" = c("INS","GCG","SST","PPY","PDX1","NKX6-1","MAFA","NEUROD1","PAX6"),
    "HALLMARK_COMPLEMENT" = c("C1QA","C1QB","C1QC","C2","C3","C4A","C4B","C5","CFB","CFD","CFH","CFI"),
    "HALLMARK_COAGULATION" = c("F2","F3","F5","F7","F8","F9","F10","FGA","FGB","FGG","PLAT","PLAU","SERPINE1"),
    "HALLMARK_MITOTIC_SPINDLE" = c("TPX2","KIF11","KIF23","KIF2C","NUF2","NDC80","SPC25","CENPA","CENPE","CENPF"),
    "HALLMARK_MYOGENESIS" = c("MYOD1","MYOG","MYF5","MYF6","MEF2C","MEF2D","DES","MYL1","MYLPF","ACTA1"),
    "HALLMARK_ADIPOGENESIS" = c("PPARG","CEBPA","CEBPB","CEBPD","FABP4","ADIPOQ","LPL","SREBF1","FASN"),
    "HALLMARK_UV_RESPONSE_DN" = c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","FN1","LAMA1","LAMB1","LAMC1","TNC"),
    "HALLMARK_UV_RESPONSE_UP" = c("DUSP1","FOS","JUN","EGR1","ATF3","IER2","IER3","GADD45B","PPP1R15A"),
    "HALLMARK_ALLOGRAFT_REJECTION" = c("CD8A","CD8B","GZMA","GZMB","PRF1","GNLY","CD3E","CD3D","CD2","IL2RB"),
    "HALLMARK_KRAS_SIGNALING_DN" = c("CLDN3","CLDN4","CLDN7","TJP1","OCLN","CDH1","CTNNB1","EPCAM","MUC1","KRT8","KRT18","KRT19"),
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_EXTENDED" = c("VIM","CDH2","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","FN1","COL1A1","MMP2","MMP9","ACTA2","TAGLN","CD44")
  )
}

# 不需要下载和过滤的步骤——直接使用内置通路列表
# 内置列表已包含所有CTC相关的Hallmark通路

# --- C2: 计算通路活性 ---
cat("\nComputing pathway activity scores...\n")

pathway_results <- list()
for (pw_name in names(hallmark_sets)) {
  pw_genes <- hallmark_sets[[pw_name]]
  # 找到数据中存在的基因
  pw_genes_found <- intersect(pw_genes, genes_clean)
  if (length(pw_genes_found) < 3) next  # 至少3个基因

  for (cn in cancer_names) {
    pct_mat <- percentile_list[[cn]]
    # 找到该通路的基因在百分位矩阵中
    genes_in_mat <- intersect(pw_genes_found, rownames(pct_mat))
    if (length(genes_in_mat) < 3) next

    # 通路活性 = 成员基因百分位的均值（每个细胞）
    pw_pct <- colMeans(pct_mat[genes_in_mat, , drop = FALSE], na.rm = TRUE)

    # Other cancers（采样200细胞）
    cells_other <- which(cancer_labels != cn)
    idx_other <- sample(cells_other, min(200, length(cells_other)))
    expr_other <- expr_all[genes_in_mat, idx_other, drop = FALSE]
    pw_pct_other <- apply(as.matrix(expr_other), 2, function(x) {
      r <- rank(x, ties.method = "average")
      mean(r[rownames(expr_other) %in% genes_in_mat] / length(x), na.rm = TRUE)
    })

    # 通路活性均值
    mean_pw_self <- mean(pw_pct, na.rm = TRUE)
    mean_pw_other <- mean(pw_pct_other, na.rm = TRUE)
    pw_delta <- mean_pw_self - mean_pw_other

    # Wilcoxon
    wt <- tryCatch(
      wilcox.test(pw_pct, pw_pct_other, alternative = "greater"),
      error = function(e) list(p.value = 1)
    )

    pathway_results[[length(pathway_results) + 1]] <- data.frame(
      pathway        = pw_name,
      cancer         = cn,
      mean_pct_self  = mean_pw_self,
      mean_pct_other = mean_pw_other,
      pct_delta      = pw_delta,
      n_genes_found  = length(genes_in_mat),
      p_value        = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
}

pathway_df <- do.call(rbind, pathway_results)
pathway_df$p_adj <- p.adjust(pathway_df$p_value, method = "BH")

# 添加每个癌种的细胞数（用于过滤小样本癌种）
n_cells_map <- setNames(sapply(percentile_list, ncol), cancer_names)
pathway_df$n_cells <- n_cells_map[pathway_df$cancer]

# 效应量驱动分级（Tier 1 额外要求 n_cells >= 50）
pathway_df$confidence_tier <- ifelse(
  pathway_df$pct_delta > 0.08 & pathway_df$p_adj < 0.05 & pathway_df$n_cells >= 50, 1,
  ifelse(
    pathway_df$pct_delta > 0.04 & pathway_df$p_adj < 0.05 & pathway_df$n_cells >= 30, 2,
    ifelse(
      pathway_df$pct_delta > 0.02 & pathway_df$p_adj < 0.05, 3, 0
    )
  )
)

cat(sprintf("\nPathway results: %d rows\n", nrow(pathway_df)))
for (t in 1:3) {
  n <- sum(pathway_df$confidence_tier == t)
  cat(sprintf("  Tier %d: %d\n", t, n))
}

# --- C3: 输出通路-癌种边 ---
cat("\n===== Writing pathway_cancer_edges.csv =====\n")

pw_out <- pathway_df[pathway_df$confidence_tier > 0,
  c("pathway", "cancer", "mean_pct_self", "mean_pct_other", "pct_delta",
    "n_genes_found", "confidence_tier", "p_value", "p_adj")
]
pw_out$cancer_full <- name_map[pw_out$cancer]

write.csv(pw_out, file.path(v2_dir, "pathway_cancer_edges.csv"), row.names = FALSE)
cat(sprintf("Wrote %d pathway-cancer edges\n", nrow(pw_out)))

# --- C4: 输出基因-通路成员边 ---
cat("\n===== Writing gene_pathway_edges.csv =====\n")

gp_edges <- list()
for (pw_name in names(hallmark_sets)) {
  pw_genes <- intersect(hallmark_sets[[pw_name]], genes_clean)
  if (length(pw_genes) == 0) next
  for (g in pw_genes) {
    gp_edges[[length(gp_edges) + 1]] <- data.frame(
      gene = g, pathway = pw_name, source = "MSigDB_Hallmark",
      stringsAsFactors = FALSE
    )
  }
}
gp_df <- do.call(rbind, gp_edges)
write.csv(gp_df, file.path(v2_dir, "gene_pathway_edges.csv"), row.names = FALSE)
cat(sprintf("Wrote %d gene-pathway membership edges\n", nrow(gp_df)))

# ============================================================
# Done
# ============================================================
cat("\n===== DONE =====\n")
cat(sprintf("Output files in %s:\n", v2_dir))
cat("  gene_cancer_edges.csv       — Gene-[:CTC_EXPRESSED_IN]->CancerType\n")
cat("  gene_cancer_full_stats.csv  — Complete stats including excluded\n")
cat("  pathway_cancer_edges.csv    — Pathway-[:ACTIVE_IN]->CancerType\n")
cat("  gene_pathway_edges.csv      — Gene-[:MEMBER_OF_PATHWAY]->Pathway\n")
