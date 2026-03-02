#!/usr/bin/env bash
set -euo pipefail

# Child-only masking test:
# - Query is one child sample
# - Reference labels are population labels from non-child samples

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CHR="${1:-22}"
THREADS="${2:-${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN || echo 1)}}}"
FLARE_MEM_GB="${FLARE_MEM_GB:-64}"
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

RAW_VCF="${RAW_VCF:-${ROOT_DIR}/chr${CHR}/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz}"
QUERY_LIST="${QUERY_LIST:-${ROOT_DIR}/data/FLARE/test/query_child.test.txt}"
REF_PANEL="${REF_PANEL:-${ROOT_DIR}/data/FLARE/flare_ref_panel.child_only.population.txt}"
FLARE_JAR="${FLARE_JAR:-${ROOT_DIR}/tools/flare/flare.jar}"

WORK_DIR="${WORK_DIR:-${ROOT_DIR}/work/FLARE/test_chr${CHR}}"
OUT_PREFIX="${OUT_PREFIX:-${ROOT_DIR}/out/FLARE/child_only_population/test_chr${CHR}}"
mkdir -p "${WORK_DIR}" "$(dirname "${OUT_PREFIX}")"

if [[ -n "${MAP_FILE:-}" ]]; then
  MAP_FILE="${MAP_FILE}"
elif [[ -f "${ROOT_DIR}/maps/plink/plink.chr${CHR}.GRCh38.map" ]]; then
  MAP_FILE="${ROOT_DIR}/maps/plink/plink.chr${CHR}.GRCh38.map"
elif [[ -f "${ROOT_DIR}/maps/plink/chr_in_chrom_field/plink.chrchr${CHR}.GRCh38.map" ]]; then
  MAP_FILE="${ROOT_DIR}/maps/plink/chr_in_chrom_field/plink.chrchr${CHR}.GRCh38.map"
elif [[ -f "${ROOT_DIR}/maps/plink.chr${CHR}.GRCh38.map" ]]; then
  MAP_FILE="${ROOT_DIR}/maps/plink.chr${CHR}.GRCh38.map"
else
  MAP_FILE=""
fi

[[ -f "${RAW_VCF}" ]] || { echo "Missing VCF: ${RAW_VCF}" >&2; exit 1; }
[[ -f "${QUERY_LIST}" ]] || { echo "Missing query list: ${QUERY_LIST}" >&2; exit 1; }
[[ -f "${REF_PANEL}" ]] || { echo "Missing FLARE panel: ${REF_PANEL}" >&2; exit 1; }
if [[ -z "${MAP_FILE}" || ! -f "${MAP_FILE}" ]]; then
  echo "Missing FLARE per-chromosome genetic map for chr${CHR}." >&2
  echo "Run: bash ${ROOT_DIR}/code/scripts/prepare_genetic_maps.sh" >&2
  echo "Or set MAP_FILE=/absolute/path/to/plink.chr${CHR}.GRCh38.map" >&2
  exit 1
fi
[[ -f "${FLARE_JAR}" ]] || { echo "Missing jar: ${FLARE_JAR}" >&2; exit 1; }

REF_SAMPLES="${WORK_DIR}/ref.samples.txt"
REF_FILT="${WORK_DIR}/ref.pass.biallelic.snps.nomiss.dedup.vcf.gz"
QUERY_FILT="${WORK_DIR}/query.pass.biallelic.snps.nomiss.dedup.vcf.gz"
REF_COMMON_RAW="${WORK_DIR}/ref.common.raw.vcf.gz"
QUERY_COMMON_RAW="${WORK_DIR}/query.common.raw.vcf.gz"
REF_COMMON="${WORK_DIR}/ref.common.vcf.gz"
QUERY_COMMON="${WORK_DIR}/query.common.vcf.gz"
FORCE_PREP="${FORCE_PREP:-0}"

awk '{print $1}' "${REF_PANEL}" > "${REF_SAMPLES}"

need_prep=0
if [[ "${FORCE_PREP}" == "1" || ! -f "${REF_COMMON}" || ! -f "${QUERY_COMMON}" ]]; then
  need_prep=1
else
  # If cached files still contain duplicate positions, rebuild automatically.
  if [[ -n "$("${BCFTOOLS_BIN}" query -f '%CHROM\t%POS\n' "${REF_COMMON}" | sort | uniq -d | head -n 1)" ]]; then
    echo "Detected duplicate positions in cached ref.common.vcf.gz; rebuilding..." >&2
    need_prep=1
  fi
  if [[ -n "$("${BCFTOOLS_BIN}" query -f '%CHROM\t%POS\n' "${QUERY_COMMON}" | sort | uniq -d | head -n 1)" ]]; then
    echo "Detected duplicate positions in cached query.common.vcf.gz; rebuilding..." >&2
    need_prep=1
  fi
fi

if [[ "${need_prep}" == "1" ]]; then
  # Restrict to PASS biallelic SNPs, deduplicate same-position SNPs, and intersect ref/query sites.
  # Avoid site-level no-missing filter here because with large reference sets it can drop all sites.
  "${BCFTOOLS_BIN}" view --threads "${THREADS}" -S "${REF_SAMPLES}" -f .,PASS -m2 -M2 -v snps -Ou "${RAW_VCF}" \
    | "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -d snps -Oz -o "${REF_FILT}"
  "${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${REF_FILT}"

  "${BCFTOOLS_BIN}" view --threads "${THREADS}" -S "${QUERY_LIST}" -f .,PASS -m2 -M2 -v snps -Ou "${RAW_VCF}" \
    | "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -d snps -Oz -o "${QUERY_FILT}"
  "${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${QUERY_FILT}"

  "${BCFTOOLS_BIN}" isec --threads "${THREADS}" -n=2 -w1 -Oz -o "${REF_COMMON_RAW}" "${REF_FILT}" "${QUERY_FILT}"
  "${BCFTOOLS_BIN}" isec --threads "${THREADS}" -n=2 -w2 -Oz -o "${QUERY_COMMON_RAW}" "${REF_FILT}" "${QUERY_FILT}"
  "${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${REF_COMMON_RAW}"
  "${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${QUERY_COMMON_RAW}"

  # One more dedup pass after intersection to guarantee unique marker positions for FLARE.
  "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -d snps -Oz -o "${REF_COMMON}" "${REF_COMMON_RAW}"
  "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -d snps -Oz -o "${QUERY_COMMON}" "${QUERY_COMMON_RAW}"
  "${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${REF_COMMON}"
  "${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${QUERY_COMMON}"
else
  echo "Using existing prepared VCFs in ${WORK_DIR} (set FORCE_PREP=1 to rebuild)."
fi

# Final safety check before FLARE to avoid Java-side duplicate marker crash.
dup_ref="$("${BCFTOOLS_BIN}" query -f '%CHROM\t%POS\n' "${REF_COMMON}" | sort | uniq -d | head -n 1 || true)"
dup_query="$("${BCFTOOLS_BIN}" query -f '%CHROM\t%POS\n' "${QUERY_COMMON}" | sort | uniq -d | head -n 1 || true)"
if [[ -n "${dup_ref}" || -n "${dup_query}" ]]; then
  echo "Duplicate marker positions remain after preprocessing." >&2
  [[ -n "${dup_ref}" ]] && echo "Example ref duplicate: ${dup_ref}" >&2
  [[ -n "${dup_query}" ]] && echo "Example query duplicate: ${dup_query}" >&2
  echo "Try rerun with FORCE_PREP=1." >&2
  exit 1
fi

# Guard against accidental empty result sets.
ref_sites="$("${BCFTOOLS_BIN}" index -n "${REF_COMMON}")"
qry_sites="$("${BCFTOOLS_BIN}" index -n "${QUERY_COMMON}")"
if [[ "${ref_sites}" -eq 0 || "${qry_sites}" -eq 0 ]]; then
  echo "Prepared FLARE VCF has zero sites (ref=${ref_sites}, query=${qry_sites})." >&2
  echo "This indicates preprocessing filters were too strict for this chromosome." >&2
  exit 1
fi

echo "Running FLARE with prepared inputs:"
echo "  ref=${REF_COMMON}"
echo "  gt=${QUERY_COMMON}"
echo "  threads=${THREADS}"
echo "  java_heap_gb=${FLARE_MEM_GB}"

java -XX:ActiveProcessorCount="${THREADS}" -Xmx"${FLARE_MEM_GB}"g -jar "${FLARE_JAR}" \
  ref="${REF_COMMON}" \
  ref-panel="${REF_PANEL}" \
  gt="${QUERY_COMMON}" \
  map="${MAP_FILE}" \
  nthreads="${THREADS}" \
  out="${OUT_PREFIX}"

echo "FLARE test completed: ${OUT_PREFIX} (threads=${THREADS})"
