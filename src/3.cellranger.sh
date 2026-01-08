#!/usr/bin/env bash
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --partition=DCU
#SBATCH --job-name=cellranger
#SBATCH --output=logs/cellranger_%A_%a.log

set -euo pipefail

echo "[INFO] SLURM job started"
echo "[INFO] Job ID        : ${SLURM_JOB_ID}"
echo "[INFO] Array Task ID : ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Hostname      : $(hostname)"
echo "[INFO] Start time    : $(date)"

# -----------------------------
# Fix working directory (absolute root)
# -----------------------------
PROJECT_ROOT="$(pwd -P)"
echo "[INFO] PROJECT_ROOT  : ${PROJECT_ROOT}"

# Ensure logs dir exists (SLURM output path uses it)
mkdir -p "${PROJECT_ROOT}/logs"

# -----------------------------
# Cell Ranger configuration (absolute paths)
# -----------------------------
CELLRANGER="/public/home/GENE_proc/rongzz/Tools/Single_Cells/cellranger-10.0.0/bin/cellranger"
TRANSCRIPTOME="/public/home/GENE_proc/rongzz/Database/Resource/Homo_sapiens/GRCh38/single_cell/cell_ranger/refdata-gex-GRCh38-2024-A"

if [[ ! -x "${CELLRANGER}" ]]; then
  echo "[ERROR] cellranger not executable or not found: ${CELLRANGER}"
  exit 1
fi
if [[ ! -d "${TRANSCRIPTOME}" ]]; then
  echo "[ERROR] transcriptome folder not found: ${TRANSCRIPTOME}"
  exit 1
fi

echo "[INFO] cellranger path : ${CELLRANGER}"
"${CELLRANGER}" --version || true
echo "[INFO] transcriptome   : ${TRANSCRIPTOME}"

# -----------------------------
# Directory configuration (absolute paths)
# -----------------------------
data_dir="${PROJECT_ROOT}/data"
fastq_dir="${data_dir}/fastq"
output_dir="${PROJECT_ROOT}/qc"
cellranger_dir="${output_dir}/cellranger"
stage_dir="${output_dir}/cellranger_fastq_staging"

mkdir -p "${cellranger_dir}" "${stage_dir}"

# -----------------------------
# Get SRR ID for this array task (same logic as your fasterq script)
# -----------------------------
csv_file="${PROJECT_ROOT}/meta/SraRunTable.csv"
if [[ ! -f "${csv_file}" ]]; then
  echo "[ERROR] CSV file not found: ${csv_file}"
  exit 1
fi

# Skip header (line 1)
csv_line=$((SLURM_ARRAY_TASK_ID + 1))
srr_id="$(sed -n "${csv_line}p" "${csv_file}" | cut -d',' -f1 | tr -d '\r')"

if [[ -z "${srr_id}" ]]; then
  echo "[ERROR] Failed to extract SRR ID from ${csv_file} (line ${csv_line})"
  exit 1
fi

sample_name="${srr_id}"

echo "=============================================="
echo "[INFO] Processing SRR accession: ${srr_id}"
echo "[INFO] Sample name (--id/--sample): ${sample_name}"

# -----------------------------
# Locate FASTQs produced by fasterq-dump (absolute paths)
# Confirmed mapping:
#   *_2.fastq.gz = barcode+UMI  -> Cell Ranger R1
#   *_3.fastq.gz = cDNA         -> Cell Ranger R2
# -----------------------------
sample_fastq_dir="${fastq_dir}/${srr_id}"

fq_r1="${sample_fastq_dir}/${srr_id}_2.fastq.gz"
fq_r2="${sample_fastq_dir}/${srr_id}_3.fastq.gz"
fq_other="${sample_fastq_dir}/${srr_id}_1.fastq.gz"

if [[ ! -f "${fq_r1}" ]]; then
  echo "[ERROR] Barcode+UMI FASTQ not found: ${fq_r1}"
  exit 1
fi
if [[ ! -f "${fq_r2}" ]]; then
  echo "[ERROR] cDNA FASTQ not found: ${fq_r2}"
  exit 1
fi

echo "[INFO] FASTQ inputs:"
echo "[INFO]   R1 (barcode+UMI) : ${fq_r1}"
echo "[INFO]   R2 (cDNA)        : ${fq_r2}"
if [[ -f "${fq_other}" ]]; then
  echo "[INFO]   Other (sample)   : ${fq_other} (ignored)"
fi

# -----------------------------
# Stage FASTQs with Cell Ranger naming (symlink; absolute)
# -----------------------------
sample_stage_dir="${stage_dir}/${sample_name}"
mkdir -p "${sample_stage_dir}"

ln -sf "${fq_r1}" "${sample_stage_dir}/${sample_name}_S1_L001_R1_001.fastq.gz"
ln -sf "${fq_r2}" "${sample_stage_dir}/${sample_name}_S1_L001_R2_001.fastq.gz"

echo "[INFO] Staging directory prepared: ${sample_stage_dir}"
ls -lh "${sample_stage_dir}"

# -----------------------------
# Run cellranger count
# -----------------------------
echo "[INFO] Running cellranger count..."

SCRATCH_BASE="${SLURM_TMPDIR:-/tmp}"
export TMPDIR="${SCRATCH_BASE}/cellranger_tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"
echo "[INFO] TMPDIR = ${TMPDIR}"

# Avoid overwrite
if [[ -d "${cellranger_dir}/${sample_name}" ]]; then
  echo "[WARN] Output directory already exists: ${cellranger_dir}/${sample_name}"
  echo "[WARN] Skipping to avoid overwrite."
  exit 0
fi

(
  cd "${cellranger_dir}"
  "${CELLRANGER}" count \
    --id="${sample_name}" \
    --create-bam=true \
    --transcriptome="${TRANSCRIPTOME}" \
    --fastqs="${sample_stage_dir}" \
    --sample="${sample_name}" \
    --chemistry=auto \
    --localcores="${SLURM_CPUS_PER_TASK}" \
    --localmem=60
)

echo "[INFO] cellranger count completed for ${sample_name}"
echo "[INFO] Output: ${cellranger_dir}/${sample_name}"
echo "[INFO] Task ${SLURM_ARRAY_TASK_ID} finished successfully at $(date)"
