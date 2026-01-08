#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: 7.seurat_sct.R <input_h5> <scdbl_tsv> <barcodes_final> <out_dir> <sample_name> <threads>")
}
input_h5 <- args[1]
scdbl_tsv <- args[2]
barcodes_file <- args[3]
out_dir <- args[4]
sample_name <- args[5]
threads <- suppressWarnings(as.integer(args[6]))

suppressPackageStartupMessages({
  library(Seurat)
  library(DropletUtils)
  library(SummarizedExperiment)
  library(Matrix)
})

sce <- DropletUtils::read10xCounts(input_h5, type = "HDF5")
counts <- SummarizedExperiment::assay(sce, "counts")
barcodes <- NULL
if ("Barcode" %in% colnames(SummarizedExperiment::colData(sce))) {
  barcodes <- SummarizedExperiment::colData(sce)$Barcode
}
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}
if (!is.null(barcodes) && length(barcodes) == ncol(counts) && !anyDuplicated(barcodes)) {
  colnames(counts) <- barcodes
}
seu <- Seurat::CreateSeuratObject(counts = counts, project = sample_name, assay = "RNA")

gene_symbols <- NULL
if ("Symbol" %in% colnames(SummarizedExperiment::rowData(sce))) {
  gene_symbols <- SummarizedExperiment::rowData(sce)$Symbol
} else if ("symbol" %in% colnames(SummarizedExperiment::rowData(sce))) {
  gene_symbols <- SummarizedExperiment::rowData(sce)$symbol
} else {
  gene_symbols <- rownames(sce)
}
if (is.null(gene_symbols)) {
  gene_symbols <- rownames(counts)
}

mito_genes <- grepl("^MT-", gene_symbols)
lib_size <- Matrix::colSums(counts)
percent_mt <- rep(0, length(lib_size))
names(percent_mt) <- colnames(counts)
if (any(mito_genes) && any(lib_size > 0)) {
  mito_sum <- Matrix::colSums(counts[mito_genes, , drop = FALSE])
  percent_mt[lib_size > 0] <- mito_sum[lib_size > 0] / lib_size[lib_size > 0] * 100
}
seu[["percent.mt"]] <- percent_mt

scdbl <- read.table(
  scdbl_tsv,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!"barcode" %in% colnames(scdbl)) {
  stop("scDblFinder TSV is missing required 'barcode' column")
}

barcodes <- colnames(seu)
scdbl_barcodes <- scdbl$barcode

if (anyDuplicated(scdbl_barcodes)) {
  stop("Duplicate barcodes found in scDblFinder TSV")
}

normalize_barcode <- function(x) {
  tolower(gsub("[-_]", "", x))
}

barcodes_norm <- normalize_barcode(barcodes)
scdbl_norm <- normalize_barcode(scdbl_barcodes)

if (anyDuplicated(scdbl_norm)) {
  stop("Ambiguous barcode match after normalization; duplicates in scDblFinder TSV")
}

match_idx <- match(barcodes_norm, scdbl_norm)
missing <- sum(is.na(match_idx))
if (missing > 0) {
  cat(sprintf("[WARN] Dropping %d cells without scDblFinder calls\n", missing))
  keep_cells <- !is.na(match_idx)
  seu <- subset(seu, cells = barcodes[keep_cells])
  barcodes <- colnames(seu)
  barcodes_norm <- barcodes_norm[keep_cells]
  if (length(barcodes) == 0) {
    stop("No cells remaining after matching scDblFinder barcodes")
  }
  match_idx <- match(barcodes_norm, scdbl_norm)
}

scdbl_meta <- scdbl[match_idx, setdiff(colnames(scdbl), "barcode"), drop = FALSE]
rownames(scdbl_meta) <- barcodes
seu <- Seurat::AddMetaData(seu, scdbl_meta)

if (!"scDblFinder.class" %in% colnames(seu@meta.data)) {
  stop("scDblFinder.class not found in metadata after merge")
}

singlets <- seu$scDblFinder.class == "singlet"
cat(sprintf("[INFO] Singlets kept=%d, removed=%d\n", sum(singlets), sum(!singlets)))
if (sum(singlets) == 0) {
  stop("No singlets remain after scDblFinder filtering")
}
seu <- subset(seu, cells = colnames(seu)[singlets])

if (nzchar(barcodes_file) && file.exists(barcodes_file)) {
  keep <- readLines(barcodes_file)
  keep <- keep[keep %in% colnames(seu)]
  if (length(keep) == 0) {
    stop("No cells remain after applying barcodes_final filter")
  }
  seu <- subset(seu, cells = keep)
}

set.seed(123)
if (!is.na(threads) && threads > 1) {
  options(future.globals.maxSize = 8 * 1024^3)
}
seu <- Seurat::SCTransform(seu, verbose = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(out_dir, paste0(sample_name, "_sct.rds"))
saveRDS(seu, out_rds)
cat(sprintf("[INFO] Saved Seurat object: %s\n", out_rds))
