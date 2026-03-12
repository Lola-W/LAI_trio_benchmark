#!/usr/bin/env bash
set -euo pipefail

# Download/source-fetch helper only. No pip usage.
# Use with conda env from code/envs/lai-tools.conda.yml

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TOOLS_DIR="${ROOT_DIR}/tools"
mkdir -p "${TOOLS_DIR}"

require_cmd() {
  local c="$1"
  if ! command -v "${c}" >/dev/null 2>&1; then
    echo "Missing required command: ${c}" >&2
    echo "Activate/install conda env first:" >&2
    echo "  conda env create -f ${ROOT_DIR}/code/envs/lai-tools.conda.yml" >&2
    echo "  conda activate lai-benchmark-tools" >&2
    exit 1
  fi
}

require_cmd git
require_cmd curl

echo "[1/2] FLARE"
mkdir -p "${TOOLS_DIR}/flare"
FLARE_JAR_PATH="${TOOLS_DIR}/flare/flare.jar"
if [[ -n "${FLARE_JAR_URL:-}" ]]; then
  echo "Using FLARE_JAR_URL=${FLARE_JAR_URL}"
  curl -L "${FLARE_JAR_URL}" -o "${FLARE_JAR_PATH}"
else
  echo "FLARE_JAR_URL not set, trying default FLARE URLs..."
  flare_ok=0
  for u in \
    "https://faculty.washington.edu/browning/flare.jar" \
    "https://faculty.washington.edu/browning/flare/flare.jar"
  do
    if curl -L -f "${u}" -o "${FLARE_JAR_PATH}"; then
      flare_ok=1
      echo "Downloaded FLARE from: ${u}"
      break
    fi
  done
  if [[ "${flare_ok}" -ne 1 ]]; then
    echo "Could not auto-download FLARE jar." >&2
    echo "Rerun with explicit URL, e.g.:" >&2
    echo "  FLARE_JAR_URL='https://.../flare.jar' bash code/scripts/download_tools.sh" >&2
    exit 1
  fi
fi
echo "Saved: ${FLARE_JAR_PATH}"

echo "[2/2] RFMix"
RFMIX_GIT_URL="${RFMIX_GIT_URL:-https://github.com/slowkoni/rfmix}"
if [[ ! -d "${TOOLS_DIR}/rfmix/.git" ]]; then
  git clone --recursive "${RFMIX_GIT_URL}" "${TOOLS_DIR}/rfmix"
else
  echo "Exists: ${TOOLS_DIR}/rfmix (skip clone)"
fi

if [[ -f "${TOOLS_DIR}/rfmix/CMakeLists.txt" ]]; then
  echo "Detected CMake-based RFMix build."
  require_cmd cmake
  mkdir -p "${TOOLS_DIR}/rfmix/build"
  cmake -S "${TOOLS_DIR}/rfmix" -B "${TOOLS_DIR}/rfmix/build" -DCMAKE_BUILD_TYPE=Release
  cmake --build "${TOOLS_DIR}/rfmix/build" -j
elif [[ -f "${TOOLS_DIR}/rfmix/configure.ac" || -f "${TOOLS_DIR}/rfmix/Makefile.am" ]]; then
  echo "Detected Autotools-based RFMix build."
  require_cmd autoreconf
  require_cmd make
  (
    cd "${TOOLS_DIR}/rfmix"
    autoreconf --force --install
    ./configure
    make -j
  )
  mkdir -p "${TOOLS_DIR}/rfmix/build"
  if [[ -x "${TOOLS_DIR}/rfmix/rfmix" ]]; then
    ln -sf ../rfmix "${TOOLS_DIR}/rfmix/build/rfmix"
  fi
else
  echo "Could not detect RFMix build system in ${TOOLS_DIR}/rfmix" >&2
  exit 1
fi

if [[ -x "${TOOLS_DIR}/rfmix/build/rfmix" ]]; then
  echo "Built binary: ${TOOLS_DIR}/rfmix/build/rfmix"
elif [[ -x "${TOOLS_DIR}/rfmix/rfmix" ]]; then
  echo "Built binary: ${TOOLS_DIR}/rfmix/rfmix"
else
  echo "RFMix build completed but binary was not found." >&2
  exit 1
fi

echo "FLARE/RFMix download/build complete."
