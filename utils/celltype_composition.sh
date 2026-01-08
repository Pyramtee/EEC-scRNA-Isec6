#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Cell type composition job started"
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

R_SCRIPT="${SCRIPT_DIR}/celltype_composition.R"
if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "[ERROR] R script not found: ${R_SCRIPT}" >&2
  exit 1
fi

out_dir="${PROJECT_ROOT}/analysis/celltype_composition"
Rscript "${R_SCRIPT}" "${PROJECT_ROOT}" "${out_dir}"

echo "[INFO] Cell type composition outputs saved: ${out_dir}"
