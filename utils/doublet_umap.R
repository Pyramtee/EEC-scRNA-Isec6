#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: doublet_umap.R <scdbl_rds> <out_dir> <sample_name> [dims]")
}
scdbl_rds <- args[1]
out_dir <- args[2]
sample_name <- args[3]
dims <- suppressWarnings(as.integer(if (length(args) >= 4) args[4] else "30"))
if (is.na(dims) || dims < 5) dims <- 30

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(ggplot2)
})

plot_font <- Sys.getenv("PLOT_FONT", "Helvetica")
plot_size <- suppressWarnings(as.numeric(Sys.getenv("PLOT_BASE_SIZE", "8")))
if (is.na(plot_size) || plot_size <= 0) plot_size <- 8
label_size <- suppressWarnings(as.numeric(Sys.getenv("PLOT_LABEL_SIZE", "")))
if (is.na(label_size) || label_size <= 0) {
  label_size <- max(2.5, plot_size * 0.35)
}
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

sce <- readRDS(scdbl_rds)
counts <- SummarizedExperiment::assay(sce, "counts")
seu <- Seurat::CreateSeuratObject(counts = counts, project = sample_name, assay = "RNA")

if ("scDblFinder.class" %in% colnames(SummarizedExperiment::colData(sce))) {
  seu$scDblFinder.class <- SummarizedExperiment::colData(sce)$scDblFinder.class
} else {
  stop("scDblFinder.class not found in scDblFinder RDS")
}

set.seed(123)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = max(50, dims), verbose = FALSE)
seu <- RunUMAP(seu, dims = seq_len(dims), verbose = FALSE)

p <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "scDblFinder.class",
  label = TRUE,
  repel = TRUE,
  label.size = label_size
) +
  ggplot2::scale_color_manual(
    values = c(singlet = "#4DBBD5FF", doublet = "#E64B35FF"),
    na.value = "grey80"
  ) +
  ggplot2::scale_fill_manual(
    values = c(singlet = "#4DBBD5FF", doublet = "#E64B35FF"),
    na.value = "grey80"
  ) +
  plot_theme

sample_dir <- file.path(out_dir, sample_name)
dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(sample_dir, paste0(sample_name, "_umap_doublet.png"))

ggsave(out_png, plot = p, width = 5, height = 4, units = "in", dpi = plot_dpi)
