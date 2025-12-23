#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Initializing conda environment..."
source /home1/rongzz/miniconda3/etc/profile.d/conda.sh
conda activate tool_sra-tools
echo "[INFO] Conda environment 'tool_sra-tools' activated."

# -----------------------------
# Directory configuration
# -----------------------------
output_dir="../data"
sra_dir="${output_dir}/sra"

echo "[INFO] Creating output directories if not exist..."
mkdir -p "${sra_dir}"

# -----------------------------
# Main loop
# -----------------------------
echo "[INFO] Start processing SraRunTable.csv"

while read -r line; do

    # Skip header
    if [[ $line == Run,* ]]; then
        echo "[INFO] Skipping header line."
        continue
    fi

    srr_id=$(echo "$line" | cut -d',' -f1)

    echo "=============================================="
    echo "[INFO] Processing SRR accession: ${srr_id}"

    # -----------------------------
    # Step 1: prefetch
    # -----------------------------
    echo "[INFO] Downloading ${srr_id} with prefetch..."
    prefetch -O "${sra_dir}" -X 30G "${srr_id}"
    echo "[INFO] prefetch completed: ${srr_id}"


done < SraRunTable.csv

echo "[INFO] All SRR accessions downloaded successfully."
