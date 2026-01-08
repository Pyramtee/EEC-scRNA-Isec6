#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "."
out_tsv <- if (length(args) >= 2 && nzchar(args[2])) args[2] else ""
variable_genes <- suppressWarnings(as.integer(Sys.getenv("VARIABLE_GENES", "3000")))
if (is.na(variable_genes) || variable_genes <= 0) {
  stop("VARIABLE_GENES must be a positive integer")
}

suppressPackageStartupMessages({
  library(DropletUtils)
  library(Matrix)
})

qc_root <- file.path(project_root, "qc", "isec")
cb_root <- file.path(project_root, "qc", "cellbender")

count_paths <- list.files(qc_root, pattern = "counts_summary\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(count_paths) == 0) {
  stop("No counts_summary.tsv found under qc/isec")
}

summary_df <- do.call(rbind, lapply(count_paths, function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE)
}))
summary_df$sample_id <- as.character(summary_df$sample_id)

get_counts <- function(sample_id) {
  h5_path <- file.path(cb_root, sample_id, "filtered_feature_bc_matrix_cellbender.h5")
  if (!file.exists(h5_path)) {
    stop(sprintf("Missing cellbender filtered H5 for %s: %s", sample_id, h5_path))
  }
  sce <- DropletUtils::read10xCounts(h5_path, type = "HDF5")
  if ("Barcode" %in% colnames(SummarizedExperiment::colData(sce))) {
    barcodes <- SummarizedExperiment::colData(sce)$Barcode
    if (length(barcodes) == ncol(sce) && !anyDuplicated(barcodes)) {
      colnames(sce) <- barcodes
    }
  }
  counts <- SummarizedExperiment::assay(sce, "counts")
  list(sce = sce, counts = counts)
}

get_gene_symbols <- function(sce) {
  if ("Symbol" %in% colnames(SummarizedExperiment::rowData(sce))) {
    return(SummarizedExperiment::rowData(sce)$Symbol)
  }
  if ("symbol" %in% colnames(SummarizedExperiment::rowData(sce))) {
    return(SummarizedExperiment::rowData(sce)$symbol)
  }
  rownames(sce)
}

rows <- lapply(summary_df$sample_id, function(sample_id) {
  barcodes_path <- file.path(qc_root, sample_id, "barcodes_final.txt")
  if (!file.exists(barcodes_path)) {
    stop(sprintf("Missing barcodes_final for %s", sample_id))
  }
  barcodes_final <- readLines(barcodes_path)
  barcodes_final <- barcodes_final[nzchar(barcodes_final)]

  data <- get_counts(sample_id)
  counts <- data$counts
  sce <- data$sce

  keep <- intersect(colnames(counts), barcodes_final)
  if (length(keep) == 0) {
    stop(sprintf("No matching barcodes for %s", sample_id))
  }
  counts <- counts[, keep, drop = FALSE]

  lib_size <- Matrix::colSums(counts)
  detected <- Matrix::colSums(counts > 0)

  gene_symbols <- get_gene_symbols(sce)
  mito_genes <- grepl("^MT-", gene_symbols)
  mito_sum <- if (any(mito_genes)) Matrix::colSums(counts[mito_genes, , drop = FALSE]) else rep(0, ncol(counts))
  percent_mito <- ifelse(lib_size > 0, mito_sum / lib_size * 100, 0)

  list(
    sample_id = sample_id,
    min_umi = min(lib_size),
    max_umi = max(lib_size),
    min_detected_genes = min(detected),
    max_detected_genes = max(detected),
    max_percent_mito = max(percent_mito),
    cells_passed_qc = ncol(counts)
  )
})

qc_df <- do.call(rbind, lapply(rows, as.data.frame))

final_df <- merge(qc_df, summary_df, by = "sample_id", all.x = TRUE)
final_df$cells_from_cellranger <- final_df$n_cellranger
final_df$cells_filtered_out <- final_df$cells_from_cellranger - final_df$cells_passed_qc
final_df$variable_genes <- variable_genes

final_df <- final_df[, c(
  "sample_id",
  "min_umi",
  "max_umi",
  "min_detected_genes",
  "max_detected_genes",
  "max_percent_mito",
  "cells_from_cellranger",
  "cells_filtered_out",
  "cells_passed_qc",
  "variable_genes"
)]

final_df <- final_df[order(final_df$sample_id), ]

if (nzchar(out_tsv)) {
  write.table(final_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(final_df, sep = "\t", quote = FALSE, row.names = FALSE)
}
