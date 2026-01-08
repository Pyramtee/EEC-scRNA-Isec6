#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: 9.singler_annot.R <input_rds> <input_h5> <out_dir> <sample_name>")
}
input_rds <- args[1]
input_h5 <- args[2]
out_dir <- args[3]
sample_name <- args[4]

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(DropletUtils)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scuttle)
  library(Matrix)
  library(ggplot2)
})

cache_dir <- Sys.getenv("GYPSUM_CACHE_DIR", "")
if (nzchar(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  gypsum::cacheDirectory(cache_dir)
  cat(sprintf("[INFO] Gypsum cache: %s\n", cache_dir))
}

seu <- readRDS(input_rds)
if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("seurat_clusters not found in metadata; run clustering first")
}

sce_raw <- DropletUtils::read10xCounts(input_h5, type = "HDF5")
counts <- Seurat::GetAssayData(seu, assay = "RNA", layer = "counts")
gene_ids <- rownames(counts)

symbols <- NULL
rd <- SummarizedExperiment::rowData(sce_raw)
if ("Symbol" %in% colnames(rd)) {
  symbols <- rd$Symbol
} else if ("symbol" %in% colnames(rd)) {
  symbols <- rd$symbol
}
if (!is.null(symbols)) {
  names(symbols) <- rownames(sce_raw)
  symbols <- symbols[gene_ids]
  symbols <- ifelse(is.na(symbols) | symbols == "", gene_ids, symbols)
  dup_n <- sum(duplicated(symbols))
  if (dup_n > 0) {
    counts <- scuttle::sumCountsAcrossFeatures(counts, ids = symbols)
    counts <- Matrix::Matrix(counts, sparse = TRUE)
    cat(sprintf("[INFO] Collapsed %d duplicated gene symbols\n", dup_n))
  } else {
    rownames(counts) <- symbols
  }
  cat("[INFO] SingleR using gene symbols from H5 metadata\n")
}

sce_test <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
sce_test <- scuttle::logNormCounts(sce_test)

clusters <- as.character(seu$seurat_clusters)
names(clusters) <- colnames(seu)

hpca <- celldex::HumanPrimaryCellAtlasData()
broad <- SingleR::SingleR(
  test = sce_test,
  ref = hpca,
  labels = hpca$label.main,
  clusters = clusters,
  assay.type.test = "logcounts"
)

broad_labels <- setNames(broad$labels, rownames(broad))
if ("pruned.labels" %in% colnames(broad)) {
  broad_pruned <- setNames(broad$pruned.labels, rownames(broad))
} else {
  broad_pruned <- broad_labels
}

immune_regex <- "(T cell|B cell|NK|Monocyte|Macrophage|Dendritic|DC|Neutrophil|Granulocyte|Plasma|Mast|Immune)"
immune_clusters <- names(broad_pruned)[grepl(immune_regex, broad_pruned, ignore.case = TRUE)]

refined_labels <- broad_labels
refined_pruned <- broad_pruned
if (length(immune_clusters) > 0) {
  cells_immune <- clusters %in% immune_clusters
  monaco <- celldex::MonacoImmuneData()
  immune <- SingleR::SingleR(
    test = sce_test[, cells_immune],
    ref = monaco,
    labels = monaco$label.fine,
    clusters = clusters[cells_immune],
    assay.type.test = "logcounts"
  )
  immune_labels <- setNames(immune$labels, rownames(immune))
  if ("pruned.labels" %in% colnames(immune)) {
    immune_pruned <- setNames(immune$pruned.labels, rownames(immune))
  } else {
    immune_pruned <- immune_labels
  }
  refined_labels[names(immune_labels)] <- immune_labels
  refined_pruned[names(immune_pruned)] <- immune_pruned
}

broad_cell <- broad_labels[clusters]
names(broad_cell) <- colnames(seu)
broad_pruned_cell <- broad_pruned[clusters]
names(broad_pruned_cell) <- colnames(seu)
refined_cell <- refined_labels[clusters]
names(refined_cell) <- colnames(seu)
refined_pruned_cell <- refined_pruned[clusters]
names(refined_pruned_cell) <- colnames(seu)

seu$SingleR.broad <- broad_cell
seu$SingleR.broad.pruned <- broad_pruned_cell
seu$SingleR.refined <- refined_cell
seu$SingleR.refined.pruned <- refined_pruned_cell

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(out_dir, paste0(sample_name, "_singler.rds"))
saveRDS(seu, out_rds)

cluster_table <- data.frame(
  cluster = names(broad_labels),
  broad = unname(broad_labels),
  broad_pruned = unname(broad_pruned),
  refined = unname(refined_labels),
  refined_pruned = unname(refined_pruned),
  stringsAsFactors = FALSE
)
out_tsv <- file.path(out_dir, paste0(sample_name, "_singler_clusters.tsv"))
write.table(cluster_table, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

plot_font <- Sys.getenv("PLOT_FONT", "Helvetica")
plot_size <- suppressWarnings(as.numeric(Sys.getenv("PLOT_BASE_SIZE", "8")))
if (is.na(plot_size) || plot_size <= 0) plot_size <- 8
label_size <- suppressWarnings(as.numeric(Sys.getenv("PLOT_LABEL_SIZE", "")))
if (is.na(label_size) || label_size <= 0) {
  label_size <- max(2.5, plot_size * 0.35)
}
plot_dpi <- suppressWarnings(as.numeric(Sys.getenv("PLOT_DPI", "300")))
if (is.na(plot_dpi) || plot_dpi <= 0) plot_dpi <- 300

palette_npg <- c(
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
  group.by = "SingleR.refined.pruned",
  label = TRUE,
  repel = TRUE,
  label.size = label_size
) +
  ggplot2::scale_color_manual(
    values = grDevices::colorRampPalette(palette_npg)(
      max(1, length(unique(seu$SingleR.refined.pruned[!is.na(seu$SingleR.refined.pruned)])))
    ),
    na.value = "grey80"
  ) +
  ggplot2::scale_fill_manual(
    values = grDevices::colorRampPalette(palette_npg)(
      max(1, length(unique(seu$SingleR.refined.pruned[!is.na(seu$SingleR.refined.pruned)])))
    ),
    na.value = "grey80"
  ) +
  ggplot2::theme_classic(base_size = plot_size, base_family = plot_font) +
  ggplot2::theme(
    legend.position = "right",
    axis.title = element_text(size = plot_size),
    axis.text = element_text(size = plot_size),
    legend.title = element_text(size = plot_size),
    legend.text = element_text(size = plot_size * 0.9)
  )

out_png <- file.path(out_dir, paste0(sample_name, "_umap_singler_refined.png"))
ggplot2::ggsave(out_png, p, width = 5, height = 4, units = "in", dpi = plot_dpi)
