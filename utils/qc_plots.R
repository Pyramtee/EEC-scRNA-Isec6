#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: qc_plots.R <input_h5> <qc_barcodes> <out_dir> <sample_name>")
}
input_h5 <- args[1]
qc_barcodes <- args[2]
out_dir <- args[3]
sample_name <- args[4]

suppressPackageStartupMessages({
  library(DropletUtils)
  library(scuttle)
  library(SummarizedExperiment)
  library(Matrix)
  library(ggplot2)
})

plot_font <- Sys.getenv("PLOT_FONT", "Helvetica")
plot_size <- suppressWarnings(as.numeric(Sys.getenv("PLOT_BASE_SIZE", "8")))
if (is.na(plot_size) || plot_size <= 0) plot_size <- 8
plot_dpi <- suppressWarnings(as.numeric(Sys.getenv("PLOT_DPI", "300")))
if (is.na(plot_dpi) || plot_dpi <= 0) plot_dpi <- 300

plot_theme <- ggplot2::theme_classic(base_size = plot_size, base_family = plot_font) +
  ggplot2::theme(
    legend.position = "right",
    axis.title = element_text(size = plot_size),
    axis.text = element_text(size = plot_size),
    legend.title = element_text(size = plot_size),
    legend.text = element_text(size = plot_size * 0.9)
  )

sce <- DropletUtils::read10xCounts(input_h5, type = "HDF5")
if ("Barcode" %in% colnames(SummarizedExperiment::colData(sce))) {
  barcodes <- SummarizedExperiment::colData(sce)$Barcode
  if (length(barcodes) == ncol(sce) && !anyDuplicated(barcodes)) {
    colnames(sce) <- barcodes
  }
}

symbols <- NULL
rd <- SummarizedExperiment::rowData(sce)
if ("Symbol" %in% colnames(rd)) {
  symbols <- rd$Symbol
} else if ("symbol" %in% colnames(rd)) {
  symbols <- rd$symbol
}
if (is.null(symbols)) {
  symbols <- rownames(sce)
}

is_mito <- grepl("^MT-", symbols)
qc <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = is_mito))

qc_pass <- readLines(qc_barcodes)
qc_pass <- qc_pass[nzchar(qc_pass)]

df <- data.frame(
  barcode = colnames(sce),
  nCount_RNA = qc$sum,
  nFeature_RNA = qc$detected,
  percent_mt = qc$subsets_Mito_percent,
  qc_status = ifelse(colnames(sce) %in% qc_pass, "pass", "fail"),
  stringsAsFactors = FALSE
)

df$qc_status <- factor(df$qc_status, levels = c("pass", "fail"))

min_umi <- 500
min_genes <- 200
max_umi_cap <- 80000
max_genes_cap <- 11000

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
if (is.na(thr_lib_q995) || !is.finite(thr_lib_q995)) thr_lib_q995 <- Inf
if (is.na(thr_feat_q995) || !is.finite(thr_feat_q995)) thr_feat_q995 <- Inf
thr_lib_high <- min(thr_lib_high_mad, thr_lib_q995, max_umi_cap)
thr_feat_high <- min(thr_feat_high_mad, thr_feat_q995, max_genes_cap)

long_df <- data.frame(
  qc_status = rep(df$qc_status, times = 3),
  metric = rep(c("nCount_RNA", "nFeature_RNA", "percent_mt"), each = nrow(df)),
  value = c(df$nCount_RNA, df$nFeature_RNA, df$percent_mt),
  stringsAsFactors = FALSE
)

p_violin <- ggplot(long_df, aes(x = qc_status, y = value, fill = qc_status)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c(pass = "#4DBBD5FF", fail = "#E64B35FF")) +
  labs(title = paste0(sample_name, " QC metrics"), x = "QC status", y = "Value") +
  plot_theme

p_scatter <- ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA, color = percent_mt, alpha = qc_status)) +
  geom_point(size = 0.6) +
  scale_alpha_manual(values = c(pass = 0.6, fail = 0.2)) +
  scale_color_gradient(low = "grey80", high = "#E64B35FF") +
  geom_vline(xintercept = c(min_umi, thr_lib_high), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = c(min_genes, thr_feat_high), linetype = "dashed", color = "grey40") +
  labs(
    title = paste0(sample_name, " nCount vs nFeature"),
    x = "nCount_RNA",
    y = "nFeature_RNA",
    color = "percent.mt",
    alpha = "QC status"
  ) +
  plot_theme

sample_dir <- file.path(out_dir, sample_name)
dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

out_violin_png <- file.path(sample_dir, paste0(sample_name, "_qc_violin_metrics.png"))
out_scatter_png <- file.path(sample_dir, paste0(sample_name, "_qc_scatter_ncount_nfeature.png"))

ggsave(out_violin_png, plot = p_violin, width = 7, height = 4, units = "in", dpi = plot_dpi)
ggsave(out_scatter_png, plot = p_scatter, width = 5, height = 4, units = "in", dpi = plot_dpi)
