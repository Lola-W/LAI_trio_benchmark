#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

# Chromosomes 1–22, X, Y
CHRS=({1..22} X Y)

mkdir -p "${BASE_DIR}"

# Download README once
wget -c -P "${BASE_DIR}" \
  "${BASE_URL}/README_1kGP_phased_panel_110722.pdf"

for chr in "${CHRS[@]}"; do
    OUTDIR="${BASE_DIR}/chr${chr}"
    mkdir -p "${OUTDIR}"

    FILE="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

    echo "Downloading chr${chr}..."

    wget -c -P "${OUTDIR}" "${BASE_URL}/${FILE}"
    wget -c -P "${OUTDIR}" "${BASE_URL}/${FILE}.tbi"
done

echo "All downloads completed."
