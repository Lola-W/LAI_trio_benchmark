#!/usr/bin/env bash
set -euo pipefail

# Child-only masking test:
# - Query is one child sample
# - Reference labels are population labels from non-child samples

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CHR="${1:-22}"
THREADS="${2:-${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN || echo 1)}}}"
GNOMIX_MEM_GB="${GNOMIX_MEM_GB:-112}"
DEFAULT_LAI_ENV_PREFIX="${ROOT_DIR}/.env/lai-benchmark-tools"

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
PYTHON_BIN="${PYTHON_BIN:-python3}"
command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1 || {
  echo "Missing bcftools. Set LAI_ENV_PREFIX or BCFTOOLS_BIN." >&2
  exit 1
}
command -v "${TABIX_BIN}" >/dev/null 2>&1 || {
  echo "Missing tabix. Set LAI_ENV_PREFIX or TABIX_BIN." >&2
  exit 1
}
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || {
  echo "Missing python. Set LAI_ENV_PREFIX or PYTHON_BIN." >&2
  exit 1
}

RAW_VCF="${RAW_VCF:-${ROOT_DIR}/chr${CHR}/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz}"
QUERY_LIST="${QUERY_LIST:-${ROOT_DIR}/data/Gnomix/test/query_child.test.txt}"
SAMPLE_MAP="${SAMPLE_MAP:-${ROOT_DIR}/data/Gnomix/gnomix_sample_map.child_only.population.tsv}"
GNOMIX_SCRIPT="${GNOMIX_SCRIPT:-${ROOT_DIR}/tools/gnomix/gnomix.py}"
GNOMIX_CONFIG="${GNOMIX_CONFIG:-${ROOT_DIR}/tools/gnomix/config.yaml}"
GNOMIX_MODEL_INFERENCE="${GNOMIX_MODEL_INFERENCE:-fast}"

WORK_DIR="${WORK_DIR:-${ROOT_DIR}/work/Gnomix/test_chr${CHR}}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/out/Gnomix/child_only_population/test_chr${CHR}}"
LOG_DIR="${LOG_DIR:-${ROOT_DIR}/out/Gnomix/logs}"
mkdir -p "${WORK_DIR}" "${OUT_DIR}" "${LOG_DIR}"

RUN_TS="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/run_gnomix_chr${CHR}_${RUN_TS}.log"
exec > >(tee -a "${LOG_FILE}") 2>&1
echo "Log file: ${LOG_FILE}"
echo "Start: $(date)"
echo "Requested memory target (scheduler): ${GNOMIX_MEM_GB}G"

if [[ -n "${GNOMIX_MAP:-}" ]]; then
  GNOMIX_MAP="${GNOMIX_MAP}"
elif [[ -f "${ROOT_DIR}/maps/gnomix_map_chr${CHR}.tsv" ]]; then
  GNOMIX_MAP="${ROOT_DIR}/maps/gnomix_map_chr${CHR}.tsv"
elif [[ -f "${ROOT_DIR}/maps/plink/plink.chr${CHR}.GRCh38.map" ]]; then
  GNOMIX_MAP="${WORK_DIR}/gnomix_map_chr${CHR}.tsv"
  awk 'NF>=4{c=$1; if(c !~ /^chr/) c="chr" c; print c"\t"$4"\t"$3}' "${ROOT_DIR}/maps/plink/plink.chr${CHR}.GRCh38.map" > "${GNOMIX_MAP}"
elif [[ -f "${ROOT_DIR}/maps/plink/chr_in_chrom_field/plink.chrchr${CHR}.GRCh38.map" ]]; then
  GNOMIX_MAP="${WORK_DIR}/gnomix_map_chr${CHR}.tsv"
  awk 'NF>=4{c=$1; if(c !~ /^chr/) c="chr" c; print c"\t"$4"\t"$3}' "${ROOT_DIR}/maps/plink/chr_in_chrom_field/plink.chrchr${CHR}.GRCh38.map" > "${GNOMIX_MAP}"
elif [[ -f "${ROOT_DIR}/maps/plink.chr${CHR}.GRCh38.map" ]]; then
  GNOMIX_MAP="${WORK_DIR}/gnomix_map_chr${CHR}.tsv"
  awk 'NF>=4{c=$1; if(c !~ /^chr/) c="chr" c; print c"\t"$4"\t"$3}' "${ROOT_DIR}/maps/plink.chr${CHR}.GRCh38.map" > "${GNOMIX_MAP}"
else
  GNOMIX_MAP=""
fi

[[ -f "${RAW_VCF}" ]] || { echo "Missing VCF: ${RAW_VCF}" >&2; exit 1; }
[[ -f "${QUERY_LIST}" ]] || { echo "Missing query list: ${QUERY_LIST}" >&2; exit 1; }
[[ -f "${SAMPLE_MAP}" ]] || { echo "Missing sample map: ${SAMPLE_MAP}" >&2; exit 1; }
if [[ -z "${GNOMIX_MAP}" || ! -f "${GNOMIX_MAP}" ]]; then
  echo "Missing Gnomix genetic map for chr${CHR}." >&2
  echo "Run: bash ${ROOT_DIR}/code/scripts/prepare_genetic_maps.sh" >&2
  echo "Or set GNOMIX_MAP=/absolute/path/to/gnomix_map_chr${CHR}.tsv" >&2
  exit 1
fi
[[ -f "${GNOMIX_SCRIPT}" ]] || { echo "Missing gnomix.py: ${GNOMIX_SCRIPT}" >&2; exit 1; }
[[ -f "${GNOMIX_CONFIG}" ]] || { echo "Missing config.yaml: ${GNOMIX_CONFIG}" >&2; exit 1; }

# Keep source index fresh to avoid stale-index warnings from htslib.
if [[ ! -f "${RAW_VCF}.tbi" || "${RAW_VCF}" -nt "${RAW_VCF}.tbi" ]]; then
  echo "Refreshing index for source VCF: ${RAW_VCF}"
  "${TABIX_BIN}" -f -p vcf "${RAW_VCF}"
fi

# Auto-detect chromosome label style to avoid Gnomix reading all variants.
CONTIG_IDS="$("${BCFTOOLS_BIN}" view -h "${RAW_VCF}" | awk -F'[=,]' '/^##contig=<ID=/{print $3}')"
if printf '%s\n' "${CONTIG_IDS}" | grep -qx "chr${CHR}"; then
  GNOMIX_CHR="chr${CHR}"
elif printf '%s\n' "${CONTIG_IDS}" | grep -qx "${CHR}"; then
  GNOMIX_CHR="${CHR}"
else
  GNOMIX_CHR="chr${CHR}"
fi

REF_SAMPLES="${WORK_DIR}/ref.samples.txt"
cut -f1 "${SAMPLE_MAP}" > "${REF_SAMPLES}"

QUERY_VCF="${WORK_DIR}/query_child.vcf.gz"
REF_VCF="${WORK_DIR}/reference_child_only_population.vcf.gz"
QUERY_RAW="${WORK_DIR}/query_child.raw.vcf.gz"
REF_RAW="${WORK_DIR}/reference_child_only_population.raw.vcf.gz"

# Restrict to the requested chromosome and biallelic SNPs for a stable, bounded benchmark input size.
"${BCFTOOLS_BIN}" view --threads "${THREADS}" -r "${GNOMIX_CHR}" -S "${QUERY_LIST}" -f .,PASS -m2 -M2 -v snps -Ou "${RAW_VCF}" \
  | "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -d snps -Oz -o "${QUERY_RAW}"
"${BCFTOOLS_BIN}" view --threads "${THREADS}" -r "${GNOMIX_CHR}" -S "${REF_SAMPLES}" -f .,PASS -m2 -M2 -v snps -Ou "${RAW_VCF}" \
  | "${BCFTOOLS_BIN}" norm --threads "${THREADS}" -d snps -Oz -o "${REF_RAW}"
"${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${QUERY_RAW}"
"${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${REF_RAW}"

# Keep only shared sites between query and reference.
"${BCFTOOLS_BIN}" isec --threads "${THREADS}" -n=2 -w2 -Oz -o "${QUERY_VCF}" "${REF_RAW}" "${QUERY_RAW}"
"${BCFTOOLS_BIN}" isec --threads "${THREADS}" -n=2 -w1 -Oz -o "${REF_VCF}" "${REF_RAW}" "${QUERY_RAW}"
"${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${QUERY_VCF}"
"${TABIX_BIN}" -@ "${THREADS}" -f -p vcf "${REF_VCF}"

q_sites="$("${BCFTOOLS_BIN}" index -n "${QUERY_VCF}")"
r_sites="$("${BCFTOOLS_BIN}" index -n "${REF_VCF}")"
echo "Prepared Gnomix inputs: chrom=${GNOMIX_CHR} ref_sites=${r_sites} query_sites=${q_sites}"
if [[ "${q_sites}" -eq 0 || "${r_sites}" -eq 0 ]]; then
  echo "No common SNPs after preprocessing for ${GNOMIX_CHR}." >&2
  exit 1
fi

# Use a per-run config copy so Gnomix honors requested thread count.
THREAD_CONFIG="${WORK_DIR}/gnomix.config.threads.yaml"
awk -v n="${THREADS}" -v mode="${GNOMIX_MODEL_INFERENCE}" '
  /^[[:space:]]*model:[[:space:]]*$/ { in_model=1; print; next }
  in_model && /^[^[:space:]]/ { in_model=0 }
  {
    if (!done_cores && $0 ~ /^[[:space:]]*n_cores:[[:space:]]*/) {
      sub(/n_cores:[[:space:]]*.*/, "n_cores: " n)
      done_cores=1
    }
    if (in_model && !done_mode && $0 ~ /^[[:space:]]*inference:[[:space:]]*$/) {
      sub(/inference:[[:space:]]*.*/, "  inference: " mode)
      done_mode=1
    }
    print
  }
' "${GNOMIX_CONFIG}" > "${THREAD_CONFIG}"

# Constrain native libraries often used by ML stack.
export OMP_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

set +e
echo "Launching Gnomix (training mode from scratch, this can be long)."
echo "  threads=${THREADS}"
echo "  mem_gb_target=${GNOMIX_MEM_GB}"
"${PYTHON_BIN}" "${GNOMIX_SCRIPT}" \
  "${QUERY_VCF}" \
  "${OUT_DIR}" \
  "${GNOMIX_CHR}" \
  False \
  "${GNOMIX_MAP}" \
  "${REF_VCF}" \
  "${SAMPLE_MAP}" \
  "${THREAD_CONFIG}"
status=$?
set -e

if [[ "${status}" -ne 0 ]]; then
  if [[ "${status}" -eq 137 ]]; then
    echo "Gnomix was killed (exit 137), likely by OOM." >&2
    echo "Try: request more job memory, reduce THREADS, or lower simulation.r_admixed in config." >&2
    echo "Official alternative: use pre-trained mode for quick inference tests." >&2
  fi
  exit "${status}"
fi

if [[ ! -f "${OUT_DIR}/query_results.msp" || ! -f "${OUT_DIR}/query_results.fb" ]]; then
  echo "Gnomix exited successfully but expected outputs are missing: ${OUT_DIR}/query_results.msp / query_results.fb" >&2
  exit 1
fi

echo "Gnomix test completed: ${OUT_DIR} (threads=${THREADS}, model_inference=${GNOMIX_MODEL_INFERENCE})"
echo "End: $(date)"
