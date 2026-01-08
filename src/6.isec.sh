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

R_SCRIPT="${SCRIPT_DIR}/6.isec.R"
if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "[ERROR] R script not found: ${R_SCRIPT}"
  exit 1
fi

# -----------------------------
# Directory configuration
# -----------------------------
qc_dir="${PROJECT_ROOT}/qc"
cellranger_dir="${qc_dir}/cellranger"
cellbender_dir="${qc_dir}/cellbender"
scdbl_dir="${qc_dir}/scdblfinder"
out_root="${qc_dir}/isec"

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

for sample_name in "${samples[@]}"; do
  echo "=============================================="
  echo "[INFO] Processing SRR accession: ${sample_name}"

  cr_barcodes="${cellranger_dir}/${sample_name}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  cb_barcodes="${cellbender_dir}/${sample_name}/${sample_name}_cellbender_cell_barcodes.csv"
  scdbl_tsv="${scdbl_dir}/${sample_name}/${sample_name}_scdblfinder.tsv"
  qc_barcodes="${scdbl_dir}/${sample_name}/${sample_name}_qc_pass_barcodes.txt"

  if [[ ! -f "${cr_barcodes}" ]]; then
    echo "[ERROR] Cell Ranger barcodes not found: ${cr_barcodes}"
    exit 1
  fi
  if [[ ! -f "${cb_barcodes}" ]]; then
    echo "[ERROR] CellBender barcodes not found: ${cb_barcodes}"
    exit 1
  fi
  if [[ ! -f "${scdbl_tsv}" ]]; then
    echo "[ERROR] scDblFinder TSV not found: ${scdbl_tsv}"
    exit 1
  fi
  if [[ ! -f "${qc_barcodes}" ]]; then
    echo "[ERROR] QC barcodes not found: ${qc_barcodes}"
    exit 1
  fi

  sample_out_dir="${out_root}/${sample_name}"
  out_barcodes="${sample_out_dir}/barcodes_final.txt"
  out_summary="${sample_out_dir}/counts_summary.tsv"

  if [[ -f "${out_barcodes}" && "${FORCE:-0}" != "1" ]]; then
    echo "[WARN] Output already exists: ${out_barcodes}"
    echo "[WARN] Skipping to avoid overwrite."
    continue
  fi

  echo "[INFO] Running Isec intersection..."
  Rscript "${R_SCRIPT}" "${cr_barcodes}" "${cb_barcodes}" "${scdbl_tsv}" "${qc_barcodes}" "${sample_out_dir}" "${sample_name}"
  echo "[INFO] Isec completed for ${sample_name}"
  echo "[INFO] Output barcodes : ${out_barcodes}"
  echo "[INFO] Output summary  : ${out_summary}"

done

echo "[INFO] All samples processed."
