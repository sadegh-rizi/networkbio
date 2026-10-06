#!/usr/bin/env bash
# 05_massive_proteomics.sh — Ionescu proteomics + secretome, MassIVE MSV000095283.
#
# On 2026-09-29 the MassIVE page for MSV000095283 showed "MassIVE Private Dataset",
# so anonymous download is expected to fail until the authors release it.
# If you obtain a reviewer login from the authors, export
#   MASSIVE_USER=MSV000095283   MASSIVE_PASS=<password>
# before running. The script first checks whether the dataset is public.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
ACC=MSV000095283
out="${DATA_ROOT}/Ionescu/MassIVE/${ACC}"; mkdir -p "${out}"

page="$(curl -fsL "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=${ACC}" || true)"
if grep -q "Private Dataset" <<< "${page}" && [[ -z "${MASSIVE_USER:-}" ]]; then
  log "${ACC} is still private. Set MASSIVE_USER/MASSIVE_PASS or re-check later."; exit 1
fi

auth=(); [[ -n "${MASSIVE_USER:-}" ]] && auth=(--user="${MASSIVE_USER}" --password="${MASSIVE_PASS:-}")
# Public MassIVE datasets live under ftp://massive-ftp.ucsd.edu/vNN/<ACC>/ ; find which vNN.
root=""
for v in v01 v02 v03 v04 v05 v06 v07 v08 v09 v10 v11 v12; do
  if curl -fs --list-only ${MASSIVE_USER:+--user "${MASSIVE_USER}:${MASSIVE_PASS:-}"} "ftp://massive-ftp.ucsd.edu/${v}/${ACC}/" >/dev/null 2>&1; then
    root="ftp://massive-ftp.ucsd.edu/${v}/${ACC}/"; break
  fi
done
if [[ -z "${root}" ]]; then
  # older layout / private reviewer access
  root="ftp://massive-ftp.ucsd.edu/${ACC}/"
fi
log "listing ${root}"
curl -fs --list-only ${MASSIVE_USER:+--user "${MASSIVE_USER}:${MASSIVE_PASS:-}"} "${root}" || { log "cannot list ${root}"; exit 1; }
# Download everything except raw spectra first (ccms_peak, search results, metadata);
# pass --raw to also mirror the raw/ folder.
excl=(--reject-regex '/raw/'); [[ "${1:-}" == "--raw" ]] && excl=()
wget -c -r -np -nH --cut-dirs=2 "${auth[@]}" "${excl[@]}" -P "${out}" "${root}"
