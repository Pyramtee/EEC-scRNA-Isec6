#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: 6.isec.R <cellranger_barcodes> <cellbender_barcodes> <scdbl_tsv> [qc_barcodes] <out_dir> <sample_name>")
}
cellranger_barcodes <- args[1]
cellbender_barcodes <- args[2]
scdbl_tsv <- args[3]
has_qc <- length(args) >= 6
qc_barcodes <- if (has_qc) args[4] else ""
out_dir <- if (has_qc) args[5] else args[4]
sample_name <- if (has_qc) args[6] else args[5]

read_lines <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, "rt")
    on.exit(close(con), add = TRUE)
    return(readLines(con))
  }
  readLines(path)
}

cr <- read_lines(cellranger_barcodes)
cb <- read_lines(cellbender_barcodes)
scdbl <- read.table(scdbl_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

if (!"barcode" %in% colnames(scdbl)) {
  stop("scDblFinder TSV missing 'barcode' column")
}
class_col <- "scDblFinder.class"
if (!class_col %in% colnames(scdbl)) {
  stop("scDblFinder TSV missing 'scDblFinder.class' column")
}

cr <- unique(cr[nzchar(cr)])
cb <- unique(cb[nzchar(cb)])
scdbl_barcodes <- unique(scdbl$barcode[scdbl[[class_col]] == "singlet"])

cr_set <- unique(cr)
cb_set <- unique(cb)
scdbl_set <- unique(scdbl_barcodes)

qc_set <- NULL
if (nzchar(qc_barcodes)) {
  if (!file.exists(qc_barcodes)) {
    stop(sprintf("QC barcodes file not found: %s", qc_barcodes))
  }
  qc <- read_lines(qc_barcodes)
  qc_set <- unique(qc[nzchar(qc)])
}

isec <- cb_set[cb_set %in% cr_set & cb_set %in% scdbl_set]
if (!is.null(qc_set)) {
  isec <- isec[isec %in% qc_set]
}

if (length(isec) == 0) {
  stop("No barcodes left after intersection")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_barcodes <- file.path(out_dir, "barcodes_final.txt")
writeLines(isec, out_barcodes)

summary_df <- data.frame(
  sample_id = sample_name,
  n_cellranger = length(cr_set),
  n_cellbender = length(cb_set),
  n_scdblfinder_singlet = length(scdbl_set),
  n_qc_passed = if (is.null(qc_set)) NA_integer_ else length(qc_set),
  n_intersection = length(isec),
  stringsAsFactors = FALSE
)
out_summary <- file.path(out_dir, "counts_summary.tsv")
write.table(summary_df, out_summary, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("[INFO] Isec counts: cellranger=%d cellbender=%d scdbl_singlet=%d intersection=%d\n",
            length(cr_set), length(cb_set), length(scdbl_set), length(isec)))
