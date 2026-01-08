#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: 5.scdblfinder.R <input_h5> <out_dir> <sample_name> <threads>")
}
input_h5 <- args[1]
out_dir <- args[2]
sample_name <- args[3]
threads <- suppressWarnings(as.integer(args[4]))

suppressPackageStartupMessages({
  library(DropletUtils)
  library(scDblFinder)
  library(BiocParallel)
  library(SingleCellExperiment)
  library(scuttle)
})

sce <- DropletUtils::read10xCounts(input_h5, type = "HDF5")
if ("Barcode" %in% colnames(SummarizedExperiment::colData(sce))) {
  barcodes <- SummarizedExperiment::colData(sce)$Barcode
  if (length(barcodes) == ncol(sce) && !anyDuplicated(barcodes)) {
    colnames(sce) <- barcodes
  }
}

gene_symbols <- NULL
if ("Symbol" %in% colnames(SummarizedExperiment::rowData(sce))) {
  gene_symbols <- SummarizedExperiment::rowData(sce)$Symbol
} else if ("symbol" %in% colnames(SummarizedExperiment::rowData(sce))) {
  gene_symbols <- SummarizedExperiment::rowData(sce)$symbol
} else {
  gene_symbols <- rownames(sce)
}
if (is.null(gene_symbols)) {
  gene_symbols <- rownames(sce)
}

is_mito <- grepl("^MT-", gene_symbols)
qc <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = is_mito))
min_umi <- 500
min_genes <- 200
max_umi_cap <- 80000
max_genes_cap <- 11000
low_lib_hard <- qc$sum < min_umi
low_feat_hard <- qc$detected < min_genes
max_mito <- 15
if (any(is_mito)) {
  high_mito_fixed <- qc$subsets_Mito_percent > max_mito
} else {
  high_mito_fixed <- rep(FALSE, ncol(sce))
}

mad_high_log10 <- function(x) {
  logx <- log10(pmax(x, 1))
  med <- median(logx, na.rm = TRUE)
  madv <- stats::mad(logx, na.rm = TRUE)
  thr_log <- med + 3 * madv
  thr <- 10^thr_log
  if (!is.finite(thr)) {
    return(Inf)
  }
  thr
}

thr_lib_high_mad <- mad_high_log10(qc$sum)
thr_feat_high_mad <- mad_high_log10(qc$detected)
thr_lib_q995 <- as.numeric(stats::quantile(qc$sum, probs = 0.995, names = FALSE, na.rm = TRUE))
thr_feat_q995 <- as.numeric(stats::quantile(qc$detected, probs = 0.995, names = FALSE, na.rm = TRUE))
if (is.na(thr_lib_q995) || !is.finite(thr_lib_q995)) {
  thr_lib_q995 <- Inf
}
if (is.na(thr_feat_q995) || !is.finite(thr_feat_q995)) {
  thr_feat_q995 <- Inf
}
thr_lib_high <- min(thr_lib_high_mad, thr_lib_q995, max_umi_cap)
thr_feat_high <- min(thr_feat_high_mad, thr_feat_q995, max_genes_cap)
high_lib <- qc$sum > thr_lib_high
high_feat <- qc$detected > thr_feat_high

drop <- low_lib_hard | high_lib | low_feat_hard | high_feat | high_mito_fixed

cat(sprintf("[INFO] QC outlier filtering removed=%d\n", sum(drop)))
cat(sprintf("[INFO] QC thresholds: lib_low_hard=%s lib_high_mad_log10=%s lib_q995=%s lib_cap=%s lib_high_used=%s feat_low_hard=%s feat_high_mad_log10=%s feat_q995=%s feat_cap=%s feat_high_used=%s mito_fixed=%s\n",
            min_umi, thr_lib_high_mad, thr_lib_q995, max_umi_cap, thr_lib_high,
            min_genes, thr_feat_high_mad, thr_feat_q995, max_genes_cap, thr_feat_high,
            max_mito))
if (all(drop)) {
  stop("All cells filtered out by MAD outlier QC")
}
qc_passed <- colnames(sce)[!drop]
out_qc_barcodes <- file.path(out_dir, paste0(sample_name, "_qc_pass_barcodes.txt"))
writeLines(qc_passed, out_qc_barcodes)
cat(sprintf("[INFO] QC-passed barcodes saved: %s\n", out_qc_barcodes))
sce <- sce[, !drop]

bp <- BiocParallel::SerialParam()
bp_mode <- tolower(Sys.getenv("SCDBLFINDER_BPPARAM", "serial"))
if (bp_mode %in% c("multicore", "multi", "mc") && !is.na(threads) && threads > 1) {
  bp <- BiocParallel::MulticoreParam(workers = threads)
}
BiocParallel::register(bp, default = TRUE)
cat(sprintf("[INFO] BPPARAM: %s\n", class(bp)[1]))

set.seed(123)
sce <- scDblFinder(sce, BPPARAM = bp)

out_rds <- file.path(out_dir, paste0(sample_name, "_scdblfinder.rds"))
saveRDS(sce, out_rds)

df <- as.data.frame(SummarizedExperiment::colData(sce))
if ("Barcode" %in% colnames(df)) {
  df$barcode <- df$Barcode
  df$Barcode <- NULL
} else {
  df$barcode <- colnames(sce)
}
preferred <- c("barcode", "scDblFinder.class", "scDblFinder.score", "scDblFinder.weighted")
keep <- preferred[preferred %in% colnames(df)]
if (length(keep) == 0) {
  keep <- colnames(df)
}

out_tsv <- file.path(out_dir, paste0(sample_name, "_scdblfinder.tsv"))
write.table(df[, keep, drop = FALSE], out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
