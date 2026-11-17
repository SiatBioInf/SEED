#!/usr/bin/env Rscript

# Shared configuration for SEED database analysis scripts.

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_file <- sub("^--file=", "", file_arg[1])
    if (script_file != "-") {
      return(dirname(normalizePath(script_file)))
    }
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  return(getwd())
}

SCRIPT_DIR <- get_script_dir()
SCRIPT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_ROOT, ".."), mustWork = TRUE)

PATHS <- list(
  root = PROJECT_ROOT,
  code = SCRIPT_ROOT,
  config = SCRIPT_DIR,
  gene_sets_csv = file.path(
    PROJECT_ROOT,
    "data",
    "Supplementary Table S3. Curated gene sets and source references for the 31 functional states.csv"
  ),

  # Main analysis inputs. RDS is preferred for the R/Seurat workflow.
  sc_rds = file.path(PROJECT_ROOT, "data", "SEED_scRNA.rds"),
  bulk_rds = file.path(PROJECT_ROOT, "data", "SEED_bulkRNA.rds"),
  pbmc_rds = file.path(PROJECT_ROOT, "data", "Non_CTC_10xPBMC_n28098.rds"),

  # Exchange/release inputs. Convert H5AD files once before analysis.
  sc_h5ad = file.path(PROJECT_ROOT, "data", "SEED_scRNA.h5ad"),
  bulk_h5ad = file.path(PROJECT_ROOT, "data", "SEED_bulkRNA.h5ad"),
  pbmc_h5ad = file.path(PROJECT_ROOT, "data", "Non_CTC_10xPBMC_n28098.h5ad"),
  gene_annotation_gtf = Sys.getenv(
    "SEED_GTF",
    unset = file.path(PROJECT_ROOT, "data", "gencode.v44.annotation.gtf.gz")
  ),

  output = file.path(PROJECT_ROOT, "outputs"),
  confidence_output = file.path(PROJECT_ROOT, "outputs", "confidence_validation"),
  sc_state_output = file.path(PROJECT_ROOT, "outputs", "scRNA_state_analysis"),
  sc_cancer_output = file.path(PROJECT_ROOT, "outputs", "scRNA_cancer_vs_others"),
  bulk_cancer_output = file.path(PROJECT_ROOT, "outputs", "bulkRNA_cancer_vs_others"),
  release_output = file.path(PROJECT_ROOT, "outputs", "release_tables")
)

dir.create(PATHS$output, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$confidence_output, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$sc_state_output, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$sc_cancer_output, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$bulk_cancer_output, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$release_output, showWarnings = FALSE, recursive = TRUE)

required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Please install required packages first: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

STATE_THEMES <- data.frame(
  theme = c(
    rep("Epithelial-mesenchymal plasticity", 7),
    rep("Immune and inflammatory activity", 5),
    rep("Dormancy, reactivation and stress", 6),
    rep("Proliferation and genome maintenance", 8),
    rep("Metabolism and protein homeostasis", 5)
  ),
  score = c(
    "Epithelial_Score", "Mesenchymal_Score", "EMT_Score", "EMT_I_Score", "EMT_II_Score", "EMT_III_Score", "EMT_IV_Score",
    "Immune_Score", "IFN_Response_Score", "IFN_MHC_I_Score", "IFN_MHC_II_Score", "Inflammation_Score",
    "Dormancy_Score", "Reactivation_Score", "Hypoxia_Score", "Cell_Stress_Score", "Stress_InVitro_Score", "UPR_Score",
    "CellCycle_G2M_Score", "CellCycle_G1S_Score", "CellCycle_Histone_Score", "DNA_Repair_Score",
    "Apoptosis_Score", "Senescence_Score", "MYC_Score", "Chromatin_Score",
    "Energy_Metabolism_Score", "Metal_Response_Score", "Proteasome_Score", "ProteinFolding_Score", "Translation_Score"
  ),
  stringsAsFactors = FALSE
)

CANCER_COLORS <- c(
  BC = "#A6CEE3", NSCLC = "#1F78B4", SCLC = "#B2DF8A", CRC = "#33A02C",
  GC = "#FB9A99", HCC = "#E31A1C", PDAC = "#FDBF6F", PC = "#FF7F00",
  HNSC = "#CAB2D6", NPC = "#6A3D9A", OV = "#FFFF99", MELA = "#B15928"
)

plot_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 13, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 11),
    legend.text = ggplot2::element_text(size = 9),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(18, 10, 10, 10)
  )

save_plot <- function(plot, base_path, width = 8, height = 6, dpi = 300) {
  p <- plot + plot_theme
  ggplot2::ggsave(paste0(base_path, ".jpg"), plot = p, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paste0(base_path, ".pdf"), plot = p, width = width, height = height)
  invisible(p)
}

load_seurat_object <- function(rds_path, h5ad_path = NULL) {
  if (!is.null(rds_path) && file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  if (!is.null(h5ad_path) && file.exists(h5ad_path)) {
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      stop("RDS not found and SeuratDisk is required to read h5ad: ", h5ad_path, call. = FALSE)
    }
    h5seurat_path <- sub("\\.h5ad$", ".h5seurat", h5ad_path, ignore.case = TRUE)
    if (!file.exists(h5seurat_path)) {
      SeuratDisk::Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
    }
    return(SeuratDisk::LoadH5Seurat(h5seurat_path))
  }
  stop("No valid Seurat input found. Checked RDS: ", rds_path, " and h5ad: ", h5ad_path, call. = FALSE)
}

set_default_assay <- function(object, preferred = "RNA") {
  assays <- Seurat::Assays(object)
  if (length(assays) == 0) {
    stop("No assay found in Seurat object.", call. = FALSE)
  }
  assay_use <- if (preferred %in% assays) preferred else assays[1]
  Seurat::DefaultAssay(object) <- assay_use
  if (assay_use != preferred) {
    message("Using assay '", assay_use, "' because assay '", preferred, "' was not found.")
  }
  object
}

require_rds <- function(rds_path, label = "Seurat object") {
  if (!file.exists(rds_path)) {
    stop(
      label, " RDS not found: ", rds_path,
      "\nRun scripts/02_preprocessing/01_convert_h5ad_to_rds.R first, or update PATHS in scripts/01_config/01_setup_paths.R.",
      call. = FALSE
    )
  }
  readRDS(rds_path)
}

clean_genes <- function(x) {
  x <- unique(toupper(trimws(as.character(x))))
  x <- x[!is.na(x) & x != "" & grepl("[A-Z]", x)]
  x <- x[!grepl("^MP[0-9]+\\.", x)]
  x
}

message("SEED project root: ", PROJECT_ROOT)
