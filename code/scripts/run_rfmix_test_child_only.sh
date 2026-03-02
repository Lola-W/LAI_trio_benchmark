#!/usr/bin/env bash
set -euo pipefail

# Child-only masking test:
# - Query is one child sample
# - Reference labels are population labels from non-child samples

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CHR="${1:-22}"
THREADS="${2:-${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN || echo 1)}}}"
RFMIX_THREADS="${RFMIX_THREADS:-${THREADS}}"
RFMIX_MEM_GB="${RFMIX_MEM_GB:-112}"
DEFAULT_LAI_ENV_PREFIX="/tscc/nfs/home/jiweng/ps-gleesonlab5/user/jiweng/conda-envs/lai-benchmark-tools"

LAI_ENV_PREFIX="${LAI_ENV_PREFIX:-${LAI_TOOLS_ENV_PREFIX:-}}"
if [[ -z "${LAI_ENV_PREFIX}" ]]; then
  if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/bcftools" ]]; then
    LAI_ENV_PREFIX="${CONDA_PREFIX}"
  elif [[ -x "${DEFAULT_LAI_ENV_PREFIX}/bin/bcftools" ]]; then
    LAI_ENV_PREFIX="${DEFAULT_LAI_ENV_PREFIX}"
  fi
fi
if [[ -n "${LAI_ENV_PREFIX}" ]]; then
  export PATH="${LAI_ENV_PREFIX}/bin:${PATH}"
fi
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
TABIX_BIN="${TABIX_BIN:-tabix}"
command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1 || {
  echo "Missing bcftools. Set LAI_ENV_PREFIX or BCFTOOLS_BIN." >&2
  exit 1
}
command -v "${TABIX_BIN}" >/dev/null 2>&1 || {
  echo "Missing tabix. Set LAI_ENV_PREFIX or TABIX_BIN." >&2
  exit 1
}

# Standalone defaults:
# - Use all child-masked reference samples (no population downsampling)
# - Use all biallelic SNPs passing FILTER (no extra AF thinning)
# - Use RFMix2 built-in defaults unless user explicitly overrides.
# Snakemake can override MAX_REF_PER_POP/REF_MIN_AF for OOM-safe runs.
MAX_REF_PER_POP="${MAX_REF_PER_POP:-0}"   # 0 => keep all reference samples
REF_MIN_AF="${REF_MIN_AF:-0}"             # 0 => no AF filter in reference

RAW_VCF="${RAW_VCF:-${ROOT_DIR}/chr${CHR}/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz}"
QUERY_LIST="${QUERY_LIST:-${ROOT_DIR}/data/RFMix/test/query_child.test.txt}"
SAMPLE_MAP="${SAMPLE_MAP:-${ROOT_DIR}/data/RFMix/rfmix_sample_map.child_only.population.tsv}"

if [[ -n "${RFMIX_BIN:-}" ]]; then
  RFMIX_BIN="${RFMIX_BIN}"
elif [[ -x "${ROOT_DIR}/tools/rfmix/build/rfmix" ]]; then
  RFMIX_BIN="${ROOT_DIR}/tools/rfmix/build/rfmix"
elif [[ -x "${ROOT_DIR}/tools/rfmix/rfmix" ]]; then
  RFMIX_BIN="${ROOT_DIR}/tools/rfmix/rfmix"
else
  RFMIX_BIN="${ROOT_DIR}/tools/rfmix/build/rfmix"
fi

WORK_DIR="${WORK_DIR:-${ROOT_DIR}/work/RFMix/test_chr${CHR}}"
OUT_PREFIX="${OUT_PREFIX:-${ROOT_DIR}/out/RFMix/child_only_population/test_chr${CHR}}"
LOG_DIR="${LOG_DIR:-${ROOT_DIR}/out/RFMix/logs}"
mkdir -p "${WORK_DIR}" "$(dirname "${OUT_PREFIX}")" "${LOG_DIR}"

RUN_TS="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/run_rfmix_chr${CHR}_${RUN_TS}.log"
exec > >(tee -a "${LOG_FILE}") 2>&1
echo "Log file: ${LOG_FILE}"
echo "Start: $(date)"
echo "Requested memory target (scheduler): ${RFMIX_MEM_GB}G"

if [[ -n "${GENETIC_MAP:-}" ]]; then
  GENETIC_MAP="${GENETIC_MAP}"
elif [[ -f "${ROOT_DIR}/maps/plink_allchr.GRCh38.map.tsv" ]]; then
  GENETIC_MAP="${ROOT_DIR}/maps/plink_allchr.GRCh38.map.tsv"
elif [[ -f "${ROOT_DIR}/maps/plink/plink.chr${CHR}.GRCh38.map" ]]; then
  GENETIC_MAP="${WORK_DIR}/rfmix_genetic_map_chr${CHR}.tsv"
  awk 'NF>=4{c=$1; if(c !~ /^chr/) c="chr" c; print c"\t"$4"\t"$3}' "${ROOT_DIR}/maps/plink/plink.chr${CHR}.GRCh38.map" > "${GENETIC_MAP}"
elif [[ -f "${ROOT_DIR}/maps/plink/chr_in_chrom_field/plink.chrchr${CHR}.GRCh38.map" ]]; then
  GENETIC_MAP="${WORK_DIR}/rfmix_genetic_map_chr${CHR}.tsv"
  awk 'NF>=4{c=$1; if(c !~ /^chr/) c="chr" c; print c"\t"$4"\t"$3}' "${ROOT_DIR}/maps/plink/chr_in_chrom_field/plink.chrchr${CHR}.GRCh38.map" > "${GENETIC_MAP}"
elif [[ -f "${ROOT_DIR}/maps/plink.chr${CHR}.GRCh38.map" ]]; then
  GENETIC_MAP="${WORK_DIR}/rfmix_genetic_map_chr${CHR}.tsv"
  awk 'NF>=4{c=$1; if(c !~ /^chr/) c="chr" c; print c"\t"$4"\t"$3}' "${ROOT_DIR}/maps/plink.chr${CHR}.GRCh38.map" > "${GENETIC_MAP}"
else
  GENETIC_MAP=""
fi

[[ -f "${RAW_VCF}" ]] || { echo "Missing VCF: ${RAW_VCF}" >&2; exit 1; }
[[ -f "${QUERY_LIST}" ]] || { echo "Missing query list: ${QUERY_LIST}" >&2; exit 1; }
[[ -f "${SAMPLE_MAP}" ]] || { echo "Missing sample map: ${SAMPLE_MAP}" >&2; exit 1; }
if [[ -z "${GENETIC_MAP}" || ! -f "${GENETIC_MAP}" ]]; then
  echo "Missing RFMix genetic map." >&2
  echo "Run: bash ${ROOT_DIR}/code/scripts/prepare_genetic_maps.sh" >&2
  echo "Or set GENETIC_MAP=/absolute/path/to/plink_allchr.GRCh38.map.tsv" >&2
  exit 1
fi
[[ -x "${RFMIX_BIN}" ]] || { echo "Missing/Not executable RFMix binary: ${RFMIX_BIN}" >&2; exit 1; }

# Keep source index fresh to avoid stale-index warnings from htslib.
if [[ ! -f "${RAW_VCF}.tbi" || "${RAW_VCF}" -nt "${RAW_VCF}.tbi" ]]; then
  echo "Refreshing index for source VCF: ${RAW_VCF}"
  "${TABIX_BIN}" -f -p vcf "${RAW_VCF}"
fi

SUBSET_MAP="${WORK_DIR}/sample_map.subset.tsv"
if [[ "${MAX_REF_PER_POP}" -gt 0 ]]; then
  awk -v max="${MAX_REF_PER_POP}" '{
    if (++n[$2] <= max) print $1"\t"$2
  }' "${SAMPLE_MAP}" > "${SUBSET_MAP}"
else
  cp "${SAMPLE_MAP}" "${SUBSET_MAP}"
fi

REF_SAMPLES="${WORK_DIR}/ref.samples.txt"
cut -f1 "${SUBSET_MAP}" > "${REF_SAMPLES}"

if [[ ! -s "${REF_SAMPLES}" ]]; then
  echo "Reference sample list is empty after subsetting." >&2
  exit 1
fi

QUERY_VCF="${WORK_DIR}/query_child.vcf.gz"
REF_VCF_RAW="${WORK_DIR}/reference_child_only_population.raw.vcf.gz"
REF_VCF="${WORK_DIR}/reference_child_only_population.vcf.gz"
QUERY_VCF_RAW="${WORK_DIR}/query_child.raw.vcf.gz"

# Keep only PASS/. biallelic SNPs and deduplicate for consistent marker set.
"${BCFTOOLS_BIN}" view --threads "${THREADS}" -S "${QUERY_LIST}" -f .,PASS -m2 -M2 -v snps -Ou "${RAW_VCF}" \
  | "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -d snps -Oz -o "${QUERY_VCF_RAW}"
"${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${QUERY_VCF_RAW}"

if [[ "${REF_MIN_AF}" == "0" ]]; then
  "${BCFTOOLS_BIN}" view --threads "${THREADS}" -S "${REF_SAMPLES}" -f .,PASS -m2 -M2 -v snps -Ou "${RAW_VCF}" \
    | "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -d snps -Oz -o "${REF_VCF_RAW}"
else
  "${BCFTOOLS_BIN}" view --threads "${THREADS}" -S "${REF_SAMPLES}" -f .,PASS -m2 -M2 -v snps -q "${REF_MIN_AF}:minor" -Ou "${RAW_VCF}" \
    | "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -d snps -Oz -o "${REF_VCF_RAW}"
fi
"${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${REF_VCF_RAW}"

# Intersect sites so query/reference are aligned and smaller.
"${BCFTOOLS_BIN}" isec --threads "${THREADS}" -n=2 -w2 -Oz -o "${QUERY_VCF}" "${REF_VCF_RAW}" "${QUERY_VCF_RAW}"
"${BCFTOOLS_BIN}" isec --threads "${THREADS}" -n=2 -w1 -Oz -o "${REF_VCF}" "${REF_VCF_RAW}" "${QUERY_VCF_RAW}"
"${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${QUERY_VCF}"
"${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${REF_VCF}"

q_sites="$("${BCFTOOLS_BIN}" index -n "${QUERY_VCF}")"
r_sites="$("${BCFTOOLS_BIN}" index -n "${REF_VCF}")"
echo "Prepared RFMix inputs: ref_samples=$(wc -l < "${REF_SAMPLES}") ref_sites=${r_sites} query_sites=${q_sites}"
if [[ "${q_sites}" -eq 0 || "${r_sites}" -eq 0 ]]; then
  echo "No common SNPs after preprocessing." >&2
  exit 1
fi

RFMIX_ARGS=(
  -f "${QUERY_VCF}"
  -r "${REF_VCF}"
  -m "${SUBSET_MAP}"
  -g "${GENETIC_MAP}"
  -o "${OUT_PREFIX}"
  --n-threads="${RFMIX_THREADS}"
  --chromosome="chr${CHR}"
)

if [[ -n "${RFMIX_TREES:-}" ]]; then
  RFMIX_ARGS+=(-t "${RFMIX_TREES}")
fi
if [[ -n "${RFMIX_CRF_SPACING:-}" ]]; then
  RFMIX_ARGS+=(-c "${RFMIX_CRF_SPACING}")
fi
if [[ -n "${RFMIX_RF_WINDOW_SIZE:-}" ]]; then
  RFMIX_ARGS+=(-s "${RFMIX_RF_WINDOW_SIZE}")
fi

echo "Running RFMix command:"
echo "  ${RFMIX_BIN} ${RFMIX_ARGS[*]}"
echo "  threads=${RFMIX_THREADS}"
echo "  mem_gb_target=${RFMIX_MEM_GB}"
set +e
"${RFMIX_BIN}" "${RFMIX_ARGS[@]}"
status=$?
set -e

if [[ "${status}" -ne 0 ]]; then
  echo "RFMix exited with status ${status}." >&2
  if [[ "${status}" -eq 137 ]]; then
    echo "Likely OOM kill (exit 137)." >&2
    echo "Try increasing job memory and/or reducing THREADS/RFMIX_THREADS." >&2
  fi
  exit "${status}"
fi

if [[ ! -f "${OUT_PREFIX}.msp.tsv" || ! -f "${OUT_PREFIX}.fb.tsv" ]]; then
  echo "RFMix finished but expected outputs are missing: ${OUT_PREFIX}.msp.tsv / ${OUT_PREFIX}.fb.tsv" >&2
  exit 1
fi

echo "RFMix test completed: ${OUT_PREFIX} (prep_threads=${THREADS}, rfmix_threads=${RFMIX_THREADS}, max_ref_per_pop=${MAX_REF_PER_POP}, ref_min_af=${REF_MIN_AF})"
echo "End: $(date)"
