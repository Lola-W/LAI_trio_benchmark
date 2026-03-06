#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
META_DIR="${ROOT_DIR}/code/meta"
DATA_DIR="${ROOT_DIR}/data"
TEST_TRIO_ID="${1:-TRIO0001}"

for req in \
  "${META_DIR}/trios.complete.tsv" \
  "${META_DIR}/sample_to_population.tsv" \
  "${META_DIR}/sample_to_superpopulation.tsv" \
  "${META_DIR}/sample_to_pop_superpop.tsv" \
  "${META_DIR}/flare_ref_panel.child_only.population.txt" \
  "${META_DIR}/rfmix_sample_map.child_only.population.tsv" \
  "${META_DIR}/gnomix_sample_map.child_only.population.tsv"
do
  [[ -f "${req}" ]] || { echo "Missing required file: ${req}" >&2; exit 1; }
done

mkdir -p "${DATA_DIR}/common"
mkdir -p "${DATA_DIR}/FLARE/test"
mkdir -p "${DATA_DIR}/RFMix/test"
mkdir -p "${DATA_DIR}/Gnomix/test"

cp "${META_DIR}/sample_to_population.tsv" "${DATA_DIR}/common/"
cp "${META_DIR}/sample_to_superpopulation.tsv" "${DATA_DIR}/common/"
cp "${META_DIR}/sample_to_pop_superpop.tsv" "${DATA_DIR}/common/"

cp "${META_DIR}"/flare_ref_panel.* "${DATA_DIR}/FLARE/"
cp "${META_DIR}"/rfmix_sample_map.* "${DATA_DIR}/RFMix/"
cp "${META_DIR}"/gnomix_sample_map.* "${DATA_DIR}/Gnomix/"

trio_line="$(awk -v t="${TEST_TRIO_ID}" '$1==t{print; exit}' "${META_DIR}/trios.complete.tsv")"
[[ -n "${trio_line}" ]] || { echo "Trio not found: ${TEST_TRIO_ID}" >&2; exit 1; }

child_id="$(awk -v t="${TEST_TRIO_ID}" '$1==t{print $2; exit}' "${META_DIR}/trios.complete.tsv")"
father_id="$(awk -v t="${TEST_TRIO_ID}" '$1==t{print $3; exit}' "${META_DIR}/trios.complete.tsv")"
mother_id="$(awk -v t="${TEST_TRIO_ID}" '$1==t{print $4; exit}' "${META_DIR}/trios.complete.tsv")"

for tool in FLARE RFMix Gnomix; do
  printf "%s\t%s\t%s\t%s\n" "${TEST_TRIO_ID}" "${child_id}" "${father_id}" "${mother_id}" > "${DATA_DIR}/${tool}/test/trio.test.tsv"
  printf "%s\n" "${child_id}" > "${DATA_DIR}/${tool}/test/query_child.test.txt"
done

echo "Staged tool maps and test manifests under: ${DATA_DIR}"
echo "Test trio: ${TEST_TRIO_ID}, child=${child_id}, father=${father_id}, mother=${mother_id}"
