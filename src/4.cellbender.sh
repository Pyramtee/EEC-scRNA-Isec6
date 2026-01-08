#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
echo "[INFO] PROJECT_ROOT  : ${PROJECT_ROOT}"

# -----------------------------
# Conda activation (local)
# -----------------------------
CONDA_SH="${CONDA_SH:-/home/rong/miniforge3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-tool_cellbender}"

if [[ ! -f "${CONDA_SH}" ]]; then
  echo "[ERROR] conda.sh not found: ${CONDA_SH}"
  exit 1
fi

echo "[INFO] Activating conda environment: ${CONDA_ENV}"
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"
echo "[INFO] Conda env activated: ${CONDA_ENV}"

# -----------------------------
# Cuda Activation
# -----------------------------
source /etc/profile.d/cuda12.6.sh

# -----------------------------
# Verify cellbender
# -----------------------------
if ! command -v cellbender >/dev/null 2>&1; then
  echo "[ERROR] cellbender not found in conda env tool_sc-qc."
  exit 1
fi

echo "[INFO] cellbender path : $(command -v cellbender)"
cellbender --version || true

export MPLBACKEND=Agg

# -----------------------------
# Directory configuration (absolute paths)
# -----------------------------
output_dir="${PROJECT_ROOT}/qc"
cellranger_dir="${output_dir}/cellranger"
cellbender_dir="${output_dir}/cellbender"

mkdir -p "${cellbender_dir}"

# -----------------------------
# Get SRR IDs from CSV (or pass sample IDs as args)
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

# -----------------------------
TMP_BASE="${PROJECT_ROOT}/.tmp"
mkdir -p "${TMP_BASE}"

for sample_name in "${samples[@]}"; do
  echo "=============================================="
  echo "[INFO] Processing SRR accession: ${sample_name}"

  # -----------------------------
  # Locate Cell Ranger raw matrix and metrics
  # -----------------------------
  cr_sample_dir="${cellranger_dir}/${sample_name}"
  cr_outs_dir="${cr_sample_dir}/outs"

  if [[ ! -d "${cr_outs_dir}" ]]; then
    echo "[ERROR] Cell Ranger outs directory not found: ${cr_outs_dir}"
    echo "[ERROR] Please ensure cellranger count has completed for ${sample_name}."
    exit 1
  fi

  raw_h5=""
  if [[ -f "${cr_outs_dir}/raw_feature_bc_matrix.h5" ]]; then
    raw_h5="${cr_outs_dir}/raw_feature_bc_matrix.h5"
  else
    echo "[ERROR] Raw matrix H5 not found under: ${cr_outs_dir}"
    echo "[ERROR] Expected: raw_feature_bc_matrix.h5"
    exit 1
  fi

  metrics_csv="${cr_outs_dir}/metrics_summary.csv"
  if [[ ! -f "${metrics_csv}" ]]; then
    echo "[ERROR] metrics_summary.csv not found: ${metrics_csv}"
    exit 1
  fi

  echo "[INFO] Cell Ranger raw H5  : ${raw_h5}"
  echo "[INFO] Cell Ranger metrics : ${metrics_csv}"

  echo "[INFO] CellBender parameters:"
  echo "[INFO]   epochs                 : 150"
  echo "[INFO]   model                  : full"
  echo "[INFO]   fpr                    : 0.01"
  echo "[INFO]   cuda                   : enabled"

  # -----------------------------
  # Prepare output directory
  # -----------------------------
  sample_out_dir="${cellbender_dir}/${sample_name}"
  mkdir -p "${sample_out_dir}"

  out_base="${sample_name}_cellbender"
  out_h5="${sample_out_dir}/${out_base}.h5"

  # Avoid overwrite
  if [[ -f "${out_h5}" ]]; then
    echo "[WARN] Output file already exists: ${out_h5}"
    echo "[WARN] Skipping to avoid overwrite."
    continue
  fi

  # -----------------------------
  # TMPDIR (local)
  # -----------------------------
  export TMPDIR="${TMP_BASE}/cellbender_tmp_${sample_name}"
  mkdir -p "${TMPDIR}"
  echo "[INFO] TMPDIR = ${TMPDIR}"

  # -----------------------------
  # Run CellBender remove-background (GPU)
  # -----------------------------
  echo "[INFO] Running cellbender remove-background (GPU)..."

  cmd=(cellbender remove-background
    --input "${raw_h5}"
    --output "${out_h5}"
    --model full
    --epochs 150
    --fpr 0.01
    --cuda
  )

  echo "[INFO] Command:"
  printf '[INFO]   %q\n' "${cmd[@]}"
  echo

  "${cmd[@]}"

  echo "[INFO] cellbender completed for ${sample_name}"
  echo "[INFO] Output base: ${out_h5}"

  echo "[INFO] Listing output directory: ${sample_out_dir}"
  ls -lh "${sample_out_dir}"

  filtered_h5="${sample_out_dir}/${out_base}_filtered.h5"
  if [[ -f "${filtered_h5}" ]]; then
    ln -sf "${filtered_h5}" "${sample_out_dir}/filtered_feature_bc_matrix_cellbender.h5"
    echo "[INFO] Symlink created:"
    echo "[INFO]   ${sample_out_dir}/filtered_feature_bc_matrix_cellbender.h5 -> ${filtered_h5}"
  else
    echo "[WARN] Expected filtered output not found: ${filtered_h5}"
    echo "[WARN] Please inspect CellBender logs and output files in ${sample_out_dir}."
  fi
done

echo "[INFO] All samples processed."
