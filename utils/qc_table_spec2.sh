#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] SPEC2 QC table job started"
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

R_SCRIPT="${SCRIPT_DIR}/qc_table_spec2.R"
if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "[ERROR] R script not found: ${R_SCRIPT}" >&2
  exit 1
fi

out_tsv="${PROJECT_ROOT}/analysis/Table1_QC_SPEC2.tsv"
Rscript "${R_SCRIPT}" "${PROJECT_ROOT}" "${out_tsv}"

echo "[INFO] SPEC2 QC table saved: ${out_tsv}"
