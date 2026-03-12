#!/usr/bin/env bash
set -euo pipefail

# Build an RFMix-ready input bundle from refreshed FLARE inputs.
# Default behavior stages metadata/query files only (small + fast).
# Optionally stage FLARE-prepared masked VCFs per child via symlink or copy.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CHR="${1:-22}"
DEST_DIR="${2:-${ROOT_DIR}/data/RFMix_from_FLARE/child_only_population/chr${CHR}}"

FLARE_DIR="${FLARE_DIR:-${ROOT_DIR}/data/FLARE}"
FLARE_TEST_DIR="${FLARE_TEST_DIR:-${FLARE_DIR}/test}"
FLARE_PANEL="${FLARE_PANEL:-${FLARE_DIR}/flare_ref_panel.child_only.population.txt}"
FLARE_WORK_DIR="${FLARE_WORK_DIR:-${ROOT_DIR}/work/FLARE/chr${CHR}}"

INCLUDE_MASKED_VCF="${INCLUDE_MASKED_VCF:-0}"  # 0|1
MASKED_MODE="${MASKED_MODE:-link}"             # link|copy

if [[ ! -d "${FLARE_TEST_DIR}" ]]; then
  echo "Missing FLARE test dir: ${FLARE_TEST_DIR}" >&2
  exit 1
fi
if [[ ! -f "${FLARE_PANEL}" ]]; then
  echo "Missing FLARE panel file: ${FLARE_PANEL}" >&2
  exit 1
fi
if [[ "${INCLUDE_MASKED_VCF}" != "0" && "${INCLUDE_MASKED_VCF}" != "1" ]]; then
  echo "INCLUDE_MASKED_VCF must be 0 or 1, got: ${INCLUDE_MASKED_VCF}" >&2
  exit 1
fi
if [[ "${MASKED_MODE}" != "link" && "${MASKED_MODE}" != "copy" ]]; then
  echo "MASKED_MODE must be link or copy, got: ${MASKED_MODE}" >&2
  exit 1
fi
if [[ "${INCLUDE_MASKED_VCF}" == "1" && ! -d "${FLARE_WORK_DIR}" ]]; then
  echo "Missing FLARE work dir for masked VCF staging: ${FLARE_WORK_DIR}" >&2
  exit 1
fi

TEST_OUT="${DEST_DIR}/test"
mkdir -p "${DEST_DIR}" "${TEST_OUT}"

SAMPLE_MAP_OUT="${DEST_DIR}/rfmix_sample_map.child_only.population.tsv"
awk 'NF>=2 {print $1"\t"$2}' "${FLARE_PANEL}" > "${SAMPLE_MAP_OUT}"

q_tmp="${DEST_DIR}/.query_child.test.tmp"
t_tmp="${DEST_DIR}/.trio.test.tmp"
: > "${q_tmp}"
: > "${t_tmp}"

child_count=0
trio_count=0
while IFS= read -r child; do
  src_q="${FLARE_TEST_DIR}/${child}/query_child.test.txt"
  src_t="${FLARE_TEST_DIR}/${child}/trio.test.tsv"
  [[ -f "${src_q}" ]] || continue

  out_child="${TEST_OUT}/${child}"
  mkdir -p "${out_child}"
  cp -f "${src_q}" "${out_child}/query_child.test.txt"
  head -n 1 "${src_q}" >> "${q_tmp}"
  child_count=$((child_count + 1))

  if [[ -f "${src_t}" ]]; then
    cp -f "${src_t}" "${out_child}/trio.test.tsv"
    head -n 1 "${src_t}" >> "${t_tmp}"
    trio_count=$((trio_count + 1))
  fi
done < <(find "${FLARE_TEST_DIR}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

if [[ "${child_count}" -eq 0 ]]; then
  echo "No per-child FLARE query files found under: ${FLARE_TEST_DIR}" >&2
  rm -f "${q_tmp}" "${t_tmp}"
  exit 1
fi

awk 'NF>0 && !seen[$1]++ {print $1}' "${q_tmp}" > "${TEST_OUT}/query_child.test.txt"
if [[ -s "${t_tmp}" ]]; then
  awk 'NF>0 && !seen[$2]++ {print $0}' "${t_tmp}" > "${TEST_OUT}/trio.test.tsv"
fi
rm -f "${q_tmp}" "${t_tmp}"

if [[ -f "${ROOT_DIR}/maps/plink_allchr.GRCh38.map.tsv" ]]; then
  cp -f "${ROOT_DIR}/maps/plink_allchr.GRCh38.map.tsv" "${DEST_DIR}/plink_allchr.GRCh38.map.tsv"
fi

staged_masked=0
missing_masked=0
if [[ "${INCLUDE_MASKED_VCF}" == "1" ]]; then
  masked_out="${DEST_DIR}/masked_vcf"
  mkdir -p "${masked_out}"
  while IFS= read -r child; do
    src_dir="${FLARE_WORK_DIR}/${child}"
    src_q="${src_dir}/query.common.vcf.gz"
    src_q_tbi="${src_q}.tbi"
    src_r="${src_dir}/ref.common.vcf.gz"
    src_r_tbi="${src_r}.tbi"

    if [[ ! -f "${src_q}" || ! -f "${src_q_tbi}" || ! -f "${src_r}" || ! -f "${src_r_tbi}" ]]; then
      missing_masked=$((missing_masked + 1))
      continue
    fi

    out_child="${masked_out}/${child}"
    mkdir -p "${out_child}"
    dst_q="${out_child}/query_child.vcf.gz"
    dst_r="${out_child}/reference_child_only_population.vcf.gz"

    if [[ "${MASKED_MODE}" == "copy" ]]; then
      cp -f "${src_q}" "${dst_q}"
      cp -f "${src_q_tbi}" "${dst_q}.tbi"
      cp -f "${src_r}" "${dst_r}"
      cp -f "${src_r_tbi}" "${dst_r}.tbi"
    else
      ln -sfn "${src_q}" "${dst_q}"
      ln -sfn "${src_q_tbi}" "${dst_q}.tbi"
      ln -sfn "${src_r}" "${dst_r}"
      ln -sfn "${src_r_tbi}" "${dst_r}.tbi"
    fi
    staged_masked=$((staged_masked + 1))
  done < "${TEST_OUT}/query_child.test.txt"
fi

echo "RFMix bundle staged from FLARE:"
echo "  dest_dir: ${DEST_DIR}"
echo "  sample_map: ${SAMPLE_MAP_OUT}"
echo "  query_children: ${child_count}"
echo "  trio_rows: ${trio_count}"
if [[ "${INCLUDE_MASKED_VCF}" == "1" ]]; then
  echo "  masked_vcf_staged: ${staged_masked} (mode=${MASKED_MODE})"
  echo "  masked_vcf_missing: ${missing_masked}"
fi

echo
echo "Example command for one child:"
echo "  CHILD=\$(head -n1 \"${TEST_OUT}/query_child.test.txt\")"
echo "  RAW_VCF=\"${ROOT_DIR}/chr${CHR}/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz\" \\"
echo "  QUERY_LIST=\"${TEST_OUT}/\${CHILD}/query_child.test.txt\" \\"
echo "  SAMPLE_MAP=\"${SAMPLE_MAP_OUT}\" \\"
echo "  GENETIC_MAP=\"${DEST_DIR}/plink_allchr.GRCh38.map.tsv\" \\"
echo "  bash \"${ROOT_DIR}/code/scripts/run_rfmix_test_child_only.sh\" \"${CHR}\" 14"
