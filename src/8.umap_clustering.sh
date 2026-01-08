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

R_SCRIPT="${SCRIPT_DIR}/8.umap_clustering.R"
if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "[ERROR] R script not found: ${R_SCRIPT}"
  exit 1
fi

Rscript -e 'if (!requireNamespace("Seurat", quietly=TRUE)) stop("Seurat not installed")'
Rscript -e 'if (!requireNamespace("ggplot2", quietly=TRUE)) stop("ggplot2 not installed")'

# -----------------------------
# Directory configuration
# -----------------------------
in_root="${PROJECT_ROOT}/analysis/seurat"
out_root="${PROJECT_ROOT}/analysis/umap_clustering"

mkdir -p "${out_root}"

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

dims="${UMAP_DIMS:-30}"
resolution="${CLUSTER_RESOLUTION:-0.6}"

for sample_name in "${samples[@]}"; do
  echo "=============================================="
  echo "[INFO] Processing SRR accession: ${sample_name}"

  input_rds="${in_root}/${sample_name}/${sample_name}_sct.rds"
  if [[ ! -f "${input_rds}" ]]; then
    echo "[ERROR] Input Seurat object not found: ${input_rds}"
    exit 1
  fi

  sample_out_dir="${out_root}/${sample_name}"
  out_rds="${sample_out_dir}/${sample_name}_umap.rds"
  out_png="${sample_out_dir}/${sample_name}_umap_clusters.png"

  if [[ -f "${out_rds}" && -f "${out_png}" && "${FORCE:-0}" != "1" ]]; then
    echo "[WARN] Output already exists: ${out_rds}"
    echo "[WARN] Skipping to avoid overwrite."
    continue
  fi

  echo "[INFO] Running UMAP/clustering..."
  Rscript "${R_SCRIPT}" "${input_rds}" "${sample_out_dir}" "${sample_name}" "${dims}" "${resolution}"
done

echo "[INFO] All samples processed."
