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

R_SCRIPT="${SCRIPT_DIR}/7.seurat_sct.R"
if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "[ERROR] R script not found: ${R_SCRIPT}"
  exit 1
fi

Rscript -e 'if (!requireNamespace("Seurat", quietly=TRUE)) stop("Seurat not installed")'
Rscript -e 'if (!requireNamespace("DropletUtils", quietly=TRUE)) stop("DropletUtils not installed")'
Rscript -e 'if (!requireNamespace("SummarizedExperiment", quietly=TRUE)) stop("SummarizedExperiment not installed")'

# -----------------------------
# Directory configuration
# -----------------------------
qc_dir="${PROJECT_ROOT}/qc"
cellbender_dir="${qc_dir}/cellbender"
scdbl_dir="${qc_dir}/scdblfinder"
isec_dir="${qc_dir}/isec"
out_root="${PROJECT_ROOT}/analysis/seurat"

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

threads="${THREADS:-}"
if [[ -z "${threads}" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    threads="$(nproc)"
  else
    threads="8"
  fi
fi

for sample_name in "${samples[@]}"; do
  echo "=============================================="
  echo "[INFO] Processing SRR accession: ${sample_name}"

  input_h5="${cellbender_dir}/${sample_name}/filtered_feature_bc_matrix_cellbender.h5"
  scdbl_tsv="${scdbl_dir}/${sample_name}/${sample_name}_scdblfinder.tsv"
  barcodes_final="${isec_dir}/${sample_name}/barcodes_final.txt"

  if [[ ! -f "${input_h5}" ]]; then
    echo "[ERROR] CellBender filtered matrix not found: ${input_h5}"
    exit 1
  fi
  if [[ ! -f "${scdbl_tsv}" ]]; then
    echo "[ERROR] scDblFinder TSV not found: ${scdbl_tsv}"
    exit 1
  fi
  if [[ ! -f "${barcodes_final}" ]]; then
    echo "[ERROR] Isec barcodes_final not found: ${barcodes_final}"
    exit 1
  fi

  sample_out_dir="${out_root}/${sample_name}"
  out_rds="${sample_out_dir}/${sample_name}_sct.rds"

  if [[ -f "${out_rds}" && "${FORCE:-0}" != "1" ]]; then
    echo "[WARN] Output already exists: ${out_rds}"
    echo "[WARN] Skipping to avoid overwrite."
    continue
  fi

  echo "[INFO] Running Seurat SCTransform..."
  Rscript "${R_SCRIPT}" "${input_h5}" "${scdbl_tsv}" "${barcodes_final}" "${sample_out_dir}" "${sample_name}" "${threads}"
done

echo "[INFO] All samples processed."
