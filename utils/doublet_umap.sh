#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Doublet UMAP job started"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"

CONDA_SH="${CONDA_SH:-/home/rong/miniforge3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-tool_sc-seurat}"

if [[ ! -f "${CONDA_SH}" ]]; then
  echo "[ERROR] conda.sh not found: ${CONDA_SH}" >&2
  exit 1
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[ERROR] Rscript not found in environment: ${CONDA_ENV}" >&2
  exit 1
fi

R_SCRIPT="${SCRIPT_DIR}/doublet_umap.R"
if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "[ERROR] R script not found: ${R_SCRIPT}" >&2
  exit 1
fi

csv_file="${PROJECT_ROOT}/meta/SraRunTable.csv"
if [[ ! -f "${csv_file}" ]]; then
  echo "[ERROR] CSV file not found: ${csv_file}" >&2
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
  echo "[ERROR] No sample IDs found. Pass SRR IDs as args or fix ${csv_file}." >&2
  exit 1
fi

out_root="${PROJECT_ROOT}/analysis/doublet"
mkdir -p "${out_root}"

for sample_name in "${samples[@]}"; do
  echo "[INFO] Processing ${sample_name}"
  scdbl_rds="${PROJECT_ROOT}/qc/scdblfinder/${sample_name}/${sample_name}_scdblfinder.rds"
  if [[ ! -f "${scdbl_rds}" ]]; then
    echo "[ERROR] scDblFinder RDS not found: ${scdbl_rds}" >&2
    exit 1
  fi

  Rscript "${R_SCRIPT}" "${scdbl_rds}" "${out_root}" "${sample_name}"
  echo "[INFO] Doublet UMAP saved under ${out_root}/${sample_name}"
done

echo "[INFO] Doublet UMAP job completed"
