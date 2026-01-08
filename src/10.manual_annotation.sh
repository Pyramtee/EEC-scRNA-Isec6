#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Job started"
echo "[INFO] Hostname      : $(hostname)"
echo "[INFO] Start time    : $(date)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
echo "[INFO] PROJECT_ROOT  : ${PROJECT_ROOT}"

CONDA_SH="${CONDA_SH:-/home/rong/miniforge3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-tool_sc-seurat}"

if [[ ! -f "${CONDA_SH}" ]]; then
  echo "[ERROR] conda.sh not found: ${CONDA_SH}" >&2
  exit 1
fi

echo "[INFO] Activating conda environment: ${CONDA_ENV}"
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"
echo "[INFO] Conda env activated: ${CONDA_ENV}"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[ERROR] Rscript not found in environment: ${CONDA_ENV}" >&2
  exit 1
fi

R_SCRIPT="${SCRIPT_DIR}/10.manual_annotation.R"
if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "[ERROR] R script not found: ${R_SCRIPT}" >&2
  exit 1
fi

qc_dir="${PROJECT_ROOT}/analysis"
singler_dir="${qc_dir}/singler"
manual_root="${qc_dir}/manual_annotation"
mkdir -p "${manual_root}"

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

for sample_name in "${samples[@]}"; do
  echo "=============================================="
  echo "[INFO] Processing SRR accession: ${sample_name}"

  input_rds="${singler_dir}/${sample_name}/${sample_name}_singler.rds"
  input_h5="${PROJECT_ROOT}/qc/cellbender/${sample_name}/filtered_feature_bc_matrix_cellbender.h5"
  if [[ ! -f "${input_rds}" ]]; then
    echo "[ERROR] SingleR RDS not found: ${input_rds}" >&2
    exit 1
  fi
  if [[ ! -f "${input_h5}" ]]; then
    echo "[ERROR] CellBender H5 not found: ${input_h5}" >&2
    exit 1
  fi

  sample_out_dir="${manual_root}/${sample_name}"
  mkdir -p "${sample_out_dir}"

  out_map="${sample_out_dir}/ManualAnnotation/manual_celltype_map.tsv"
  if [[ -f "${out_map}" && "${FORCE:-0}" != "1" ]]; then
    echo "[WARN] Manual annotation outputs exist for ${sample_name}."
    echo "[WARN] Set FORCE=1 to overwrite."
    continue
  fi

  Rscript "${R_SCRIPT}" "${input_rds}" "${sample_out_dir}" "${sample_name}" "${input_h5}"
  echo "[INFO] Manual annotation outputs saved: ${sample_out_dir}"
done

echo "[INFO] All samples processed."
