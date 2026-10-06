#!/usr/bin/env bash
# Shared helpers for the Ionescu 2024 / Park 2025 download scripts.
# Source this file; do not run it directly.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${SCRIPT_DIR}/manifest.tsv"
RUNS_TSV="${SCRIPT_DIR}/sra_runs.tsv"

# DATA_ROOT: where files go. Default is ./data relative to the directory you run from.
# Run from the repository root, or export DATA_ROOT=/path/to/data first.
DATA_ROOT="${DATA_ROOT:-$(pwd)/data}"
LOG_DIR="${LOG_DIR:-$(pwd)/logs}"
mkdir -p "${DATA_ROOT}" "${LOG_DIR}"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

have() { command -v "$1" >/dev/null 2>&1; }

# fetch URL DEST [EXPECTED_BYTES]
# Resumes partial downloads, writes to DEST.part and renames when complete,
# and checks the size when EXPECTED_BYTES is given.
fetch() {
  local url="$1" dest="$2" expected="${3:-}"
  [[ "${expected}" == "-" ]] && expected=""
  mkdir -p "$(dirname "${dest}")"
  if [[ -s "${dest}" ]]; then
    if [[ -z "${expected}" ]] || [[ "$(stat -c %s "${dest}")" == "${expected}" ]]; then
      log "skip (present): ${dest#"${DATA_ROOT}"/}"; return 0
    fi
    log "size mismatch, re-downloading: ${dest#"${DATA_ROOT}"/}"; mv "${dest}" "${dest}.part"
  fi
  log "get  ${url}"
  if have wget; then
    wget -q --show-progress -c --tries=5 --waitretry=10 -O "${dest}.part" "${url}"
  else
    curl -fL --retry 5 --retry-delay 10 -C - -o "${dest}.part" "${url}"
  fi
  if [[ -n "${expected}" ]]; then
    local got; got="$(stat -c %s "${dest}.part")"
    if [[ "${got}" != "${expected}" ]]; then
      log "ERROR size ${got} != expected ${expected} for ${dest}"; return 1
    fi
  fi
  mv "${dest}.part" "${dest}"
}

# fetch_tier TIER : download every manifest row whose first column is TIER
fetch_tier() {
  local tier="$1" n=0
  while IFS=$'\t' read -r -u 3 t dest url expected; do
    [[ "${t}" == "tier" || "${t}" != "${tier}" ]] && continue
    fetch "${url}" "${DATA_ROOT}/${dest}" "${expected}"
    n=$((n+1))
  done 3< "${MANIFEST}"
  log "tier ${tier}: ${n} file(s) checked"
}

# free_gb PATH : free space in GB on the filesystem holding PATH
free_gb() { df -BG --output=avail "$1" | tail -1 | tr -dc '0-9'; }

need_space() {  # need_space GB
  local need="$1" free; free="$(free_gb "${DATA_ROOT}")"
  if (( free < need )); then
    log "ERROR: need ~${need} GB free under ${DATA_ROOT}, have ${free} GB"; exit 1
  fi
  log "disk: ${free} GB free, ~${need} GB needed"
}
