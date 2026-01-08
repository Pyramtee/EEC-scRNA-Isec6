#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: 8.umap_clustering.R <input_rds> <out_dir> <sample_name> <dims> <resolution>")
}
input_rds <- args[1]
out_dir <- args[2]
sample_name <- args[3]
dims <- suppressWarnings(as.integer(args[4]))
resolution <- suppressWarnings(as.numeric(args[5]))

if (is.na(dims) || dims < 5) dims <- 30
if (is.na(resolution) || resolution <= 0) resolution <- 0.6

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

seu <- readRDS(input_rds)
DefaultAssay(seu) <- "SCT"

set.seed(123)
seu <- RunPCA(seu, npcs = max(50, dims), verbose = FALSE)
seu <- FindNeighbors(seu, dims = seq_len(dims), verbose = FALSE)
seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
seu <- RunUMAP(seu, dims = seq_len(dims), verbose = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(out_dir, paste0(sample_name, "_umap.rds"))
saveRDS(seu, out_rds)

plot_font <- Sys.getenv("PLOT_FONT", "Helvetica")
plot_size <- suppressWarnings(as.numeric(Sys.getenv("PLOT_BASE_SIZE", "8")))
if (is.na(plot_size) || plot_size <= 0) plot_size <- 8
label_size <- suppressWarnings(as.numeric(Sys.getenv("PLOT_LABEL_SIZE", "")))
if (is.na(label_size) || label_size <= 0) {
  label_size <- max(2.5, plot_size * 0.35)
}
plot_dpi <- suppressWarnings(as.numeric(Sys.getenv("PLOT_DPI", "300")))
if (is.na(plot_dpi) || plot_dpi <= 0) plot_dpi <- 300

preset_npg <- c(
  "#E64B35FF",
  "#4DBBD5FF",
  "#00A087FF",
  "#3C5488FF",
  "#F39B7FFF",
  "#8491B4FF",
  "#91D1C2FF",
  "#DC0000FF",
  "#7E6148FF",
  "#B09C85FF"
)

p <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = label_size
) +
  ggplot2::scale_color_manual(
    values = grDevices::colorRampPalette(preset_npg)(
      length(unique(seu$seurat_clusters))
    )
  ) +
  ggplot2::scale_fill_manual(
    values = grDevices::colorRampPalette(preset_npg)(
      length(unique(seu$seurat_clusters))
    )
  ) +
  ggplot2::theme_classic(base_size = plot_size, base_family = plot_font) +
  ggplot2::theme(
    legend.position = "right",
    axis.title = element_text(size = plot_size),
    axis.text = element_text(size = plot_size),
    legend.title = element_text(size = plot_size),
    legend.text = element_text(size = plot_size * 0.9)
  )

out_png <- file.path(out_dir, paste0(sample_name, "_umap_clusters.png"))
ggplot2::ggsave(out_png, p, width = 4, height = 3.5, units = "in", dpi = plot_dpi)
