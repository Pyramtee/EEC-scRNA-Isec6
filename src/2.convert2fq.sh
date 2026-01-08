#!/usr/bin/env bash
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --partition=DCU
#SBATCH --job-name=convert2fq
#SBATCH --output=logs/convert2fq_%A_%a.log

set -euo pipefail

echo "[INFO] SLURM job started"
echo "[INFO] Job ID        : ${SLURM_JOB_ID}"
echo "[INFO] Array Task ID : ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Hostname      : $(hostname)"
echo "[INFO] Start time    : $(date)"

# -----------------------------
# Conda environment
# -----------------------------
echo "[INFO] Initializing conda environment..."
source /home1/rongzz/miniconda3/etc/profile.d/conda.sh
conda activate tool_sra-tools
echo "[INFO] Conda environment 'tool_sra-tools' activated."

# -----------------------------
# Project paths
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"

# -----------------------------
# Directory configuration
# -----------------------------
output_dir="${PROJECT_ROOT}/data"
sra_dir="${output_dir}/sra"
fastq_dir="${output_dir}/fastq"

mkdir -p "${fastq_dir}"

# -----------------------------
# Get SRR ID for this array task
# -----------------------------
csv_file="${PROJECT_ROOT}/meta/SraRunTable.csv"

# Skip header (line 1)
csv_line=$((SLURM_ARRAY_TASK_ID + 1))

srr_id=$(sed -n "${csv_line}p" "${csv_file}" | cut -d',' -f1)

if [[ -z "${srr_id}" ]]; then
    echo "[ERROR] Failed to extract SRR ID from ${csv_file} (line ${csv_line})"
    exit 1
fi

echo "=============================================="
echo "[INFO] Processing SRR accession: ${srr_id}"

# -----------------------------
# Step 1: fasterq-dump
# -----------------------------
echo "[INFO] Converting ${srr_id}.sra to FASTQ..."

sample_fastq_dir="${fastq_dir}/${srr_id}"
mkdir -p "${sample_fastq_dir}"

fasterq-dump \
    --split-files \
    --include-technical \
    --threads 16 \
    --outdir "${sample_fastq_dir}" \
    "${sra_dir}/${srr_id}/${srr_id}.sra"

echo "[INFO] FASTQ generation completed for ${srr_id}"

# -----------------------------
# Step 2: compression
# -----------------------------
echo "[INFO] Compressing FASTQ files with pigz..."

pigz -p 16 "${sample_fastq_dir}/${srr_id}_1.fastq"
pigz -p 16 "${sample_fastq_dir}/${srr_id}_2.fastq"
pigz -p 16 "${sample_fastq_dir}/${srr_id}_3.fastq"

echo "[INFO] Compression completed for ${srr_id}"

echo "[INFO] Task ${SLURM_ARRAY_TASK_ID} finished successfully at $(date)"
