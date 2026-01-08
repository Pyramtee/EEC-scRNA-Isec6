#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] QC plot job started"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"

CONDA_SH="${CONDA_SH:-/home/rong/miniforge3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-tool_scDblFinder}"

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

R_SCRIPT="${SCRIPT_DIR}/qc_plots.R"
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

out_root="${PROJECT_ROOT}/analysis/qc_plots"
mkdir -p "${out_root}"

for sample_name in "${samples[@]}"; do
  echo "[INFO] Processing ${sample_name}"
  input_h5="${PROJECT_ROOT}/qc/cellbender/${sample_name}/filtered_feature_bc_matrix_cellbender.h5"
  qc_barcodes="${PROJECT_ROOT}/qc/scdblfinder/${sample_name}/${sample_name}_qc_pass_barcodes.txt"
  if [[ ! -f "${input_h5}" ]]; then
    echo "[ERROR] CellBender H5 not found: ${input_h5}" >&2
    exit 1
  fi
  if [[ ! -f "${qc_barcodes}" ]]; then
    echo "[ERROR] QC barcodes not found: ${qc_barcodes}" >&2
    exit 1
  fi

  Rscript "${R_SCRIPT}" "${input_h5}" "${qc_barcodes}" "${out_root}" "${sample_name}"
  echo "[INFO] QC plots saved under ${out_root}/${sample_name}"
done

echo "[INFO] QC plot job completed"
