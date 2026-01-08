#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "."
out_tsv <- if (length(args) >= 2 && nzchar(args[2])) args[2] else ""

suppressPackageStartupMessages({
  library(Seurat)
})

meta_path <- file.path(project_root, "meta", "SraRunTable.csv")
if (!file.exists(meta_path)) {
  stop(sprintf("Meta file not found: %s", meta_path))
}

meta <- read.csv(meta_path, stringsAsFactors = FALSE)
if (!"Run" %in% colnames(meta)) {
  stop("SraRunTable.csv missing Run column")
}

samples <- meta$Run

parse_int <- function(x) {
  if (is.na(x) || x == "") return(NA_integer_)
  as.integer(gsub(",", "", x))
}

derive_group <- function(x) {
  if (is.na(x)) return(NA_character_)
  if (grepl("normal", x, ignore.case = TRUE)) return("Control")
  if (grepl("EEC", x, ignore.case = TRUE)) return("Disease")
  x
}

stat_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(median = NA_real_, iqr = NA_real_))
  c(median = stats::median(x), iqr = stats::IQR(x))
}

rows <- lapply(samples, function(sample_id) {
  treatment <- NA_character_
  if ("treatment" %in% colnames(meta)) {
    treatment <- meta$treatment[meta$Run == sample_id][1]
  }
  group <- derive_group(treatment)

  metrics_path <- file.path(project_root, "qc", "cellranger", sample_id, "outs", "metrics_summary.csv")
  cellranger_cells <- NA_integer_
  if (file.exists(metrics_path)) {
    metrics <- read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
    if ("Estimated Number of Cells" %in% colnames(metrics)) {
      cellranger_cells <- parse_int(metrics[["Estimated Number of Cells"]][1])
    }
  }

  qc_pass_path <- file.path(project_root, "qc", "scdblfinder", sample_id, paste0(sample_id, "_qc_pass_barcodes.txt"))
  qc_pass_cells <- if (file.exists(qc_pass_path)) {
    sum(nzchar(readLines(qc_pass_path)))
  } else {
    NA_integer_
  }

  scdbl_path <- file.path(project_root, "qc", "scdblfinder", sample_id, paste0(sample_id, "_scdblfinder.tsv"))
  singlets <- NA_integer_
  doublets <- NA_integer_
  dbl_rate <- NA_real_
  if (file.exists(scdbl_path)) {
    scdbl <- read.delim(scdbl_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    if ("scDblFinder.class" %in% colnames(scdbl)) {
      singlets <- sum(scdbl$scDblFinder.class == "singlet", na.rm = TRUE)
      doublets <- sum(scdbl$scDblFinder.class == "doublet", na.rm = TRUE)
      total <- singlets + doublets
      if (total > 0) {
        dbl_rate <- doublets / total * 100
      }
    }
  }

  cb_path <- file.path(project_root, "qc", "cellbender", sample_id, paste0(sample_id, "_cellbender_cell_barcodes.csv"))
  cellbender_cells <- if (file.exists(cb_path)) {
    sum(nzchar(readLines(cb_path)))
  } else {
    NA_integer_
  }

  isec_path <- file.path(project_root, "qc", "isec", sample_id, "barcodes_final.txt")
  isec_cells <- if (file.exists(isec_path)) {
    sum(nzchar(readLines(isec_path)))
  } else {
    NA_integer_
  }

  seurat_path <- file.path(project_root, "analysis", "seurat", sample_id, paste0(sample_id, "_sct.rds"))
  ncount_stats <- c(median = NA_real_, iqr = NA_real_)
  nfeature_stats <- c(median = NA_real_, iqr = NA_real_)
  mito_stats <- c(median = NA_real_, iqr = NA_real_)
  if (file.exists(seurat_path)) {
    seu <- readRDS(seurat_path)
    md <- seu@meta.data
    if ("nCount_RNA" %in% colnames(md)) {
      ncount_stats <- stat_median_iqr(md$nCount_RNA)
    }
    if ("nFeature_RNA" %in% colnames(md)) {
      nfeature_stats <- stat_median_iqr(md$nFeature_RNA)
    }
    if ("percent.mt" %in% colnames(md)) {
      mito_stats <- stat_median_iqr(md$percent.mt)
    }
  }

  data.frame(
    sample_id = sample_id,
    group = group,
    cellranger_cells = cellranger_cells,
    qc_pass_cells = qc_pass_cells,
    scdbl_singlets = singlets,
    scdbl_doublets = doublets,
    scdbl_doublet_rate_pct = dbl_rate,
    cellbender_cells = cellbender_cells,
    isec_cells = isec_cells,
    nCount_RNA_median = ncount_stats["median"],
    nCount_RNA_IQR = ncount_stats["iqr"],
    nFeature_RNA_median = nfeature_stats["median"],
    nFeature_RNA_IQR = nfeature_stats["iqr"],
    percent_mt_median = mito_stats["median"],
    percent_mt_IQR = mito_stats["iqr"],
    stringsAsFactors = FALSE
  )
})

final_df <- do.call(rbind, rows)
final_df <- final_df[order(final_df$sample_id), ]

if (nzchar(out_tsv)) {
  write.table(final_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(final_df, sep = "\t", quote = FALSE, row.names = FALSE)
}
