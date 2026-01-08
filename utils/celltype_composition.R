#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "."
out_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else file.path(project_root, "analysis", "celltype_composition")

suppressPackageStartupMessages({
  library(Seurat)
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

pct_labels <- function(x) sprintf("%d%%", round(x * 100))

meta_path <- file.path(project_root, "meta", "SraRunTable.csv")
if (!file.exists(meta_path)) {
  stop(sprintf("Meta file not found: %s", meta_path))
}
meta <- read.csv(meta_path, stringsAsFactors = FALSE)
if (!"Run" %in% colnames(meta)) {
  stop("SraRunTable.csv missing Run column")
}

samples <- meta$Run

get_group <- function(x) {
  if (is.na(x)) return(NA_character_)
  if (grepl("normal", x, ignore.case = TRUE)) return("Control")
  if (grepl("EEC", x, ignore.case = TRUE)) return("Disease")
  x
}

rows <- list()
for (sample_id in samples) {
  rds_path <- file.path(project_root, "analysis", "singler", sample_id, paste0(sample_id, "_singler.rds"))
  map_path <- file.path(project_root, "analysis", "manual_annotation", sample_id, "ManualAnnotation", "manual_celltype_map.tsv")
  if (!file.exists(rds_path)) {
    stop(sprintf("SingleR RDS not found: %s", rds_path))
  }
  seu <- readRDS(rds_path)
  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop(sprintf("seurat_clusters missing in %s", rds_path))
  }
  if (file.exists(map_path)) {
    map <- read.delim(map_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    map <- data.frame(cluster = sort(unique(seu$seurat_clusters)), celltype_manual = "Unknown")
  }

  cluster <- as.character(seu$seurat_clusters)
  celltype <- map$celltype_manual[match(cluster, as.character(map$cluster))]
  celltype[is.na(celltype)] <- "Unknown"

  treatment <- NA_character_
  if ("treatment" %in% colnames(meta)) {
    treatment <- meta$treatment[meta$Run == sample_id][1]
  }
  group <- get_group(treatment)

  tab <- as.data.frame(table(celltype), stringsAsFactors = FALSE)
  names(tab) <- c("celltype", "n")
  tab$sample_id <- sample_id
  tab$group <- group
  rows[[sample_id]] <- tab
}

df <- do.call(rbind, rows)

sample_totals <- tapply(df$n, df$sample_id, sum)

df$proportion <- df$n / sample_totals[df$sample_id]

df$celltype <- factor(df$celltype, levels = sort(unique(df$celltype)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_sample <- ggplot(df, aes(x = sample_id, y = proportion, fill = celltype)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = pct_labels) +
  labs(title = "Cell type composition by sample", x = "Sample", y = "Proportion", fill = "Cell type") +
  plot_theme

out_sample_png <- file.path(out_dir, "celltype_composition_by_sample.png")

ggsave(out_sample_png, plot = p_sample, width = 7, height = 4.5, units = "in", dpi = plot_dpi)

agg <- aggregate(proportion ~ group + celltype, data = df, FUN = mean)

p_group <- ggplot(agg, aes(x = celltype, y = proportion, fill = group)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = pct_labels) +
  labs(title = "Cell type proportion by group", x = "Cell type", y = "Mean proportion", fill = "Group") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out_group_png <- file.path(out_dir, "celltype_composition_by_group.png")

ggsave(out_group_png, plot = p_group, width = 7.5, height = 4.5, units = "in", dpi = plot_dpi)

write.table(df, file = file.path(out_dir, "celltype_composition_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(agg, file = file.path(out_dir, "celltype_composition_by_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
