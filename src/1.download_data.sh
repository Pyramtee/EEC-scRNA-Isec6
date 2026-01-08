#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Initializing conda environment..."
source /home1/rongzz/miniconda3/etc/profile.d/conda.sh
conda activate tool_sra-tools
echo "[INFO] Conda environment 'tool_sra-tools' activated."

# -----------------------------
# Project paths
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
csv_file="${PROJECT_ROOT}/meta/SraRunTable.csv"

# -----------------------------
# Directory configuration
# -----------------------------
output_dir="${PROJECT_ROOT}/data"
sra_dir="${output_dir}/sra"

echo "[INFO] Creating output directories if not exist..."
mkdir -p "${sra_dir}"

# -----------------------------
# Main loop
# -----------------------------
if [[ ! -f "${csv_file}" ]]; then
    echo "[ERROR] CSV file not found: ${csv_file}"
    exit 1
fi

echo "[INFO] Start processing ${csv_file}"

while read -r line; do

    # Skip header
    if [[ $line == Run,* ]]; then
        echo "[INFO] Skipping header line."
        continue
    fi

    srr_id=$(echo "$line" | cut -d',' -f1 | tr -d '\r')

    echo "=============================================="
    echo "[INFO] Processing SRR accession: ${srr_id}"

    # -----------------------------
    # Step 1: prefetch
    # -----------------------------
    echo "[INFO] Downloading ${srr_id} with prefetch..."
    prefetch -O "${sra_dir}" -X 30G "${srr_id}"
    echo "[INFO] prefetch completed: ${srr_id}"


done < "${csv_file}"

echo "[INFO] All SRR accessions downloaded successfully."
