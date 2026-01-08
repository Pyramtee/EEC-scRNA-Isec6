#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Job started"
echo "[INFO] Hostname      : $(hostname)"
echo "[INFO] Start time    : $(date)"

# -----------------------------
# Project paths
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
echo "[INFO] PROJECT_ROOT  : ${PROJECT_ROOT}"

# -----------------------------
# Conda activation
# -----------------------------
CONDA_SH="${CONDA_SH:-/home/rong/miniforge3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-tool_sc-seurat}"

if [[ ! -f "${CONDA_SH}" ]]; then
  echo "[ERROR] conda.sh not found: ${CONDA_SH}"
  exit 1
fi

echo "[INFO] Activating conda environment: ${CONDA_ENV}"
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"
echo "[INFO] Conda env activated: ${CONDA_ENV}"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[ERROR] Rscript not found in environment: ${CONDA_ENV}"
  exit 1
fi

R_SCRIPT="${SCRIPT_DIR}/9.singler_annot.R"
if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "[ERROR] R script not found: ${R_SCRIPT}"
  exit 1
fi

Rscript -e 'if (!requireNamespace("Seurat", quietly=TRUE)) stop("Seurat not installed")'
Rscript -e 'if (!requireNamespace("SingleR", quietly=TRUE)) stop("SingleR not installed")'
Rscript -e 'if (!requireNamespace("celldex", quietly=TRUE)) stop("celldex not installed")'
Rscript -e 'if (!requireNamespace("DropletUtils", quietly=TRUE)) stop("DropletUtils not installed")'
Rscript -e 'if (!requireNamespace("SummarizedExperiment", quietly=TRUE)) stop("SummarizedExperiment not installed")'
Rscript -e 'if (!requireNamespace("SingleCellExperiment", quietly=TRUE)) stop("SingleCellExperiment not installed")'
Rscript -e 'if (!requireNamespace("scuttle", quietly=TRUE)) stop("scuttle not installed")'
Rscript -e 'if (!requireNamespace("ggplot2", quietly=TRUE)) stop("ggplot2 not installed")'
Rscript -e 'if (!requireNamespace("gypsum", quietly=TRUE)) stop("gypsum not installed")'

# -----------------------------
# Directory configuration
# -----------------------------
in_root="${PROJECT_ROOT}/analysis/umap_clustering"
out_root="${PROJECT_ROOT}/analysis/singler"
cellbender_dir="${PROJECT_ROOT}/qc/cellbender"

mkdir -p "${out_root}"
export GYPSUM_CACHE_DIR="${PROJECT_ROOT}/analysis/.gypsum_cache"

# -----------------------------
# Determine samples
# -----------------------------
csv_file="${PROJECT_ROOT}/meta/SraRunTable.csv"
if [[ ! -f "${csv_file}" ]]; then
  echo "[ERROR] CSV file not found: ${csv_file}"
  exit 1
fi

samples=()

if [[ $# -gt 0 ]]; then
  samples=("$@")
else
  while read -r line; do
    if [[ $line == Run,* ]]; then
      continue
    fi
    srr_id="$(echo "${line}" | cut -d',' -f1 | tr -d '\r')"
    if [[ -n "${srr_id}" ]]; then
      samples+=("${srr_id}")
    fi
  done < "${csv_file}"
fi

if [[ ${#samples[@]} -eq 0 ]]; then
  echo "[ERROR] No sample IDs found. Pass SRR IDs as args or fix ${csv_file}."
  exit 1
fi

for sample_name in "${samples[@]}"; do
  echo "=============================================="
  echo "[INFO] Processing SRR accession: ${sample_name}"

  input_rds="${in_root}/${sample_name}/${sample_name}_umap.rds"
  input_h5="${cellbender_dir}/${sample_name}/filtered_feature_bc_matrix_cellbender.h5"

  if [[ ! -f "${input_rds}" ]]; then
    echo "[ERROR] Input Seurat object not found: ${input_rds}"
    exit 1
  fi
  if [[ ! -f "${input_h5}" ]]; then
    echo "[ERROR] CellBender filtered matrix not found: ${input_h5}"
    exit 1
  fi

  sample_out_dir="${out_root}/${sample_name}"
  out_rds="${sample_out_dir}/${sample_name}_singler.rds"
  out_tsv="${sample_out_dir}/${sample_name}_singler_clusters.tsv"
  out_png="${sample_out_dir}/${sample_name}_umap_singler_refined.png"

  if [[ -f "${out_rds}" && -f "${out_tsv}" && -f "${out_png}" && "${FORCE:-0}" != "1" ]]; then
    echo "[WARN] Output already exists: ${out_rds}"
    echo "[WARN] Skipping to avoid overwrite."
    continue
  fi

  echo "[INFO] Running SingleR annotation..."
  Rscript "${R_SCRIPT}" "${input_rds}" "${input_h5}" "${sample_out_dir}" "${sample_name}"
done

echo "[INFO] All samples processed."
