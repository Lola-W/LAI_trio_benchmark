#!/usr/bin/env bash
set -euo pipefail

# Download and prepare GRCh38 genetic maps for FLARE and RFMix.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MAPS_DIR="${ROOT_DIR}/maps"
PLINK_DIR="${MAPS_DIR}/plink"
ZIP_URL="${MAP_ZIP_URL:-https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip}"
ZIP_PATH="${MAPS_DIR}/plink.GRCh38.map.zip"

mkdir -p "${MAPS_DIR}" "${PLINK_DIR}"

if [[ ! -f "${ZIP_PATH}" ]]; then
  echo "Downloading GRCh38 PLINK genetic maps..."
  curl -L -f "${ZIP_URL}" -o "${ZIP_PATH}"
else
  echo "Using existing zip: ${ZIP_PATH}"
fi

echo "Extracting maps..."
unzip -o -q "${ZIP_PATH}" -d "${PLINK_DIR}"

# Normalize placement and naming:
# write canonical files at maps/plink/plink.chr<CHR>.GRCh38.map
for chr in $(seq 1 22) X; do
  src=""
  for cand in \
    "${PLINK_DIR}/plink.chr${chr}.GRCh38.map" \
    "${PLINK_DIR}/plink.chrchr${chr}.GRCh38.map" \
    "${PLINK_DIR}/chr_in_chrom_field/plink.chrchr${chr}.GRCh38.map" \
    "${PLINK_DIR}/no_chr_in_chrom_field/plink.chr${chr}.GRCh38.map"
  do
    if [[ -f "${cand}" ]]; then
      src="${cand}"
      break
    fi
  done

  if [[ -n "${src}" ]]; then
    cp -f "${src}" "${PLINK_DIR}/plink.chr${chr}.GRCh38.map"
  fi
done

chr_count="$(ls -1 "${PLINK_DIR}"/plink.chr[0-9]*.GRCh38.map 2>/dev/null | wc -l)"
if [[ "${chr_count}" -eq 0 ]]; then
  echo "No per-chromosome PLINK maps found after unzip." >&2
  exit 1
fi

echo "Writing RFMix all-chromosome map: ${MAPS_DIR}/plink_allchr.GRCh38.map.tsv"
: > "${MAPS_DIR}/plink_allchr.GRCh38.map.tsv"
for chr in $(seq 1 22); do
  f="${PLINK_DIR}/plink.chr${chr}.GRCh38.map"
  [[ -f "${f}" ]] || continue
  # Input PLINK: chr snp_id cM bp
  # Output RFMix: chr bp cM (ensure chr prefix)
  awk 'NF>=4{c=$1; if(c !~ /^chr/) c="chr" c; print c"\t"$4"\t"$3}' "${f}" >> "${MAPS_DIR}/plink_allchr.GRCh38.map.tsv"
done

echo "Map preparation complete."
echo "FLARE per-chr maps: ${PLINK_DIR}/plink.chr<CHR>.GRCh38.map"
echo "RFMix all-chr map: ${MAPS_DIR}/plink_allchr.GRCh38.map.tsv"
