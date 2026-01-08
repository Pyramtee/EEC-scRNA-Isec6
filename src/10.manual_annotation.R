#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 10.manual_annotation.R <input_rds> <out_dir> <sample_name> [input_h5]")
}
input_rds <- args[1]
out_dir <- args[2]
sample_name <- args[3]
input_h5 <- if (length(args) >= 4) args[4] else ""

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

seu <- readRDS(input_rds)

assay_use <- if ("SCT" %in% names(seu@assays)) "SCT" else "RNA"
DefaultAssay(seu) <- assay_use

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

if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("seurat_clusters not found in metadata; run clustering first")
}

manual_dir <- file.path(out_dir, "ManualAnnotation")
heatmap_dir <- file.path(out_dir, "HeatmapEvidence")
dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

Idents(seu) <- "seurat_clusters"

plot_only <- tolower(Sys.getenv("MANUAL_PLOT_ONLY", "0")) %in% c("1", "true", "yes")

if (nzchar(input_h5) && file.exists(input_h5) && "RNA" %in% names(seu@assays)) {
  suppressPackageStartupMessages({
    library(DropletUtils)
    library(SummarizedExperiment)
    library(scuttle)
  })
  get_counts <- function(obj, assay_name) {
    tryCatch(
      SeuratObject::GetAssayData(obj, assay = assay_name, layer = "counts"),
      error = function(e) SeuratObject::GetAssayData(obj, assay = assay_name, slot = "counts")
    )
  }
  sce_raw <- DropletUtils::read10xCounts(input_h5, type = "HDF5")
  rd <- SummarizedExperiment::rowData(sce_raw)
  symbols <- NULL
  if ("Symbol" %in% colnames(rd)) {
    symbols <- rd$Symbol
  } else if ("symbol" %in% colnames(rd)) {
    symbols <- rd$symbol
  }
  if (!is.null(symbols)) {
    names(symbols) <- rownames(sce_raw)
    gene_ids <- rownames(seu@assays$RNA)
    mapped <- symbols[gene_ids]
    mapped <- ifelse(is.na(mapped) | mapped == "", gene_ids, mapped)
    counts <- get_counts(seu, "RNA")
    counts_symbol <- scuttle::sumCountsAcrossFeatures(counts, ids = mapped)
    counts_symbol <- Matrix::Matrix(counts_symbol, sparse = TRUE)
    seu[["RNA_SYMBOL"]] <- CreateAssayObject(counts_symbol)
    DefaultAssay(seu) <- "RNA_SYMBOL"
    seu <- NormalizeData(seu, assay = "RNA_SYMBOL", verbose = FALSE)
  }
}

markers <- NULL
top_markers <- NULL
marker_path <- file.path(manual_dir, "cluster_markers.tsv")
top_marker_path <- file.path(manual_dir, "top_markers_per_cluster.tsv")

if (plot_only) {
  if (file.exists(top_marker_path)) {
    top_markers <- read.delim(top_marker_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else if (file.exists(marker_path)) {
    markers <- read.delim(marker_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    if (!"gene" %in% colnames(markers)) {
      markers$gene <- rownames(markers)
    }
    fc_col <- NULL
    if ("avg_log2FC" %in% colnames(markers)) {
      fc_col <- "avg_log2FC"
    } else if ("avg_logFC" %in% colnames(markers)) {
      fc_col <- "avg_logFC"
    } else if ("avg_log" %in% colnames(markers)) {
      fc_col <- "avg_log"
    }
    if (is.null(fc_col)) {
      stop("No logFC column found in marker table")
    }
    n_top <- suppressWarnings(as.integer(Sys.getenv("TOP_MARKERS_N", "10")))
    if (is.na(n_top) || n_top <= 0) {
      n_top <- 10
    }
    split_markers <- split(markers, markers$cluster)
    top_markers <- do.call(rbind, lapply(split_markers, function(df) {
      df <- df[order(-df[[fc_col]]), , drop = FALSE]
      head(df, n_top)
    }))
  } else {
    stop("MANUAL_PLOT_ONLY=1 requires existing marker files in ManualAnnotation/")
  }
} else {
  markers <- FindAllMarkers(
    object = seu,
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.25,
    logfc.threshold = 0.25,
    return.thresh = 0.05
  )
  if (!"gene" %in% colnames(markers)) {
    markers$gene <- rownames(markers)
  }

  fc_col <- NULL
  if ("avg_log2FC" %in% colnames(markers)) {
    fc_col <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(markers)) {
    fc_col <- "avg_logFC"
  } else if ("avg_log" %in% colnames(markers)) {
    fc_col <- "avg_log"
  }
  if (is.null(fc_col)) {
    stop("No logFC column found in marker table")
  }

  markers <- markers[order(markers$cluster, -markers[[fc_col]]), ]
  write.table(
    markers,
    file = marker_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  n_top <- suppressWarnings(as.integer(Sys.getenv("TOP_MARKERS_N", "10")))
  if (is.na(n_top) || n_top <= 0) {
    n_top <- 10
  }

  split_markers <- split(markers, markers$cluster)
  top_markers <- do.call(rbind, lapply(split_markers, function(df) {
    df <- df[order(-df[[fc_col]]), , drop = FALSE]
    head(df, n_top)
  }))
  write.table(
    top_markers,
    file = top_marker_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

if (!is.null(top_markers) && !"gene" %in% colnames(top_markers)) {
  top_markers$gene <- rownames(top_markers)
}

marker_panel <- c(
  "PTPRC",
  "CD3D",
  "CD3E",
  "TRAC",
  "NKG7",
  "GNLY",
  "GZMB",
  "PRF1",
  "MS4A1",
  "CD79A",
  "CD74",
  "MZB1",
  "JCHAIN",
  "XBP1",
  "LYZ",
  "LST1",
  "FCGR3A",
  "S100A8",
  "S100A9",
  "MS4A7",
  "FCER1A",
  "CST3",
  "PECAM1",
  "VWF",
  "KDR",
  "COL1A1",
  "COL1A2",
  "DCN",
  "LUM",
  "TAGLN",
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "MKI67",
  "TOP2A",
  "HMGB2",
  "UBE2C",
  "JUN",
  "FOS",
  "HSPA1A"
)

present_features <- function(x) {
  intersect(x, rownames(seu))
}

plot_if_present <- function(features, out_path, width, height, ncol = 2) {
  features <- present_features(features)
  if (length(features) == 0) {
    message(sprintf("[WARN] No requested features found for %s", out_path))
    return(invisible(NULL))
  }
  p <- FeaturePlot(seu, features = features, ncol = ncol) + plot_theme
  ggsave(out_path, plot = p, width = width, height = height, dpi = plot_dpi)
}

plot_if_present(c("PTPRC", "EPCAM", "PECAM1", "COL1A1"), file.path(manual_dir, "featureplot_lineage_panels.png"), 10, 8, 2)
plot_if_present(c("CD3D", "MS4A1", "LYZ", "NKG7"), file.path(manual_dir, "featureplot_immune_panels.png"), 10, 8, 2)
plot_if_present(c("MKI67", "TOP2A"), file.path(manual_dir, "featureplot_cycle_panels.png"), 10, 5, 2)

panel_short <- c(
  "PTPRC","CD3D","TRAC","NKG7","GNLY","MS4A1","CD79A","LYZ","LST1",
  "PECAM1","VWF","COL1A1","DCN","EPCAM","KRT19","MKI67","TOP2A"
)

panel_short <- present_features(panel_short)
if (length(panel_short) > 0) {
  dp <- DotPlot(seu, features = panel_short) + RotatedAxis() + plot_theme
  ggsave(file.path(manual_dir, "dotplot_marker_panels.png"), plot = dp, width = 12, height = 6, dpi = plot_dpi)
} else {
  message("[WARN] No marker panel genes found for DotPlot")
}

map_path <- file.path(manual_dir, "manual_celltype_map.tsv")
map <- NULL
if (file.exists(map_path)) {
  map <- read.delim(map_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
  if ("SingleR.refined.pruned" %in% colnames(seu@meta.data)) {
    cluster_ids <- sort(unique(seu$seurat_clusters))
    celltype <- vapply(cluster_ids, function(cl) {
      labels <- seu$SingleR.refined.pruned[seu$seurat_clusters == cl]
      labels <- labels[!is.na(labels) & labels != ""]
      if (length(labels) == 0) {
        return("Unknown")
      }
      names(sort(table(labels), decreasing = TRUE))[1]
    }, character(1))
    map <- data.frame(cluster = as.character(cluster_ids), celltype_manual = unname(celltype))
  } else {
    cluster_ids <- sort(unique(seu$seurat_clusters))
    map <- data.frame(cluster = as.character(cluster_ids), celltype_manual = "Unknown")
  }
  write.table(map, file = map_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

seu$celltype_manual <- map$celltype_manual[match(as.character(seu$seurat_clusters), map$cluster)]

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

p_umap <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "celltype_manual",
  label = TRUE,
  repel = TRUE,
  label.size = label_size
) +
  ggplot2::scale_color_manual(
    values = grDevices::colorRampPalette(palette_npg)(
      max(1, length(unique(seu$celltype_manual[!is.na(seu$celltype_manual)])))
    ),
    na.value = "grey80"
  ) +
  ggplot2::scale_fill_manual(
    values = grDevices::colorRampPalette(palette_npg)(
      max(1, length(unique(seu$celltype_manual[!is.na(seu$celltype_manual)])))
    ),
    na.value = "grey80"
  ) +
  plot_theme
ggsave(file.path(manual_dir, "umap_manual_celltype.png"), plot = p_umap, width = 9, height = 7, dpi = plot_dpi)

n_down <- suppressWarnings(as.integer(Sys.getenv("HEATMAP_DOWNSAMPLE", "200")))
if (is.na(n_down) || n_down <= 0) {
  n_down <- 200
}

top_genes <- present_features(unique(top_markers$gene))
if (length(top_genes) > 0) {
  seu_small <- subset(seu, downsample = n_down)
  seu_small <- ScaleData(seu_small, features = top_genes)
  ht1 <- DoHeatmap(seu_small, features = top_genes, group.by = "seurat_clusters", raster = TRUE) + plot_theme
  ggsave(file.path(heatmap_dir, "heatmap_top_markers_per_cluster.png"), plot = ht1, width = 12, height = 10, dpi = plot_dpi)
} else {
  message("[WARN] No top marker genes found for heatmap")
}

panel_genes <- present_features(marker_panel)
if (length(panel_genes) > 0) {
  seu_small2 <- subset(seu, downsample = n_down)
  seu_small2 <- ScaleData(seu_small2, features = panel_genes)
  ht2 <- DoHeatmap(seu_small2, features = panel_genes, group.by = "celltype_manual", raster = TRUE) + plot_theme
  ggsave(file.path(heatmap_dir, "heatmap_celltype_marker_panel.png"), plot = ht2, width = 12, height = 7, dpi = plot_dpi)
} else {
  message("[WARN] No marker panel genes found for heatmap")
}
