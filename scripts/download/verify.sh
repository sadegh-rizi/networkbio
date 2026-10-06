#!/usr/bin/env bash
# verify.sh — report which manifest files are present and whether their size matches.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
ok=0; bad=0; miss=0
while IFS=$'\t' read -r t dest _url expected; do
  [[ "${t}" == "tier" ]] && continue
  f="${DATA_ROOT}/${dest}"; [[ "${expected}" == "-" ]] && expected=""
  if [[ ! -e "${f}" ]]; then miss=$((miss+1)); printf 'MISSING  %-20s %s\n' "${t}" "${dest}"; continue; fi
  if [[ -n "${expected}" && "$(stat -c %s "${f}")" != "${expected}" ]]; then bad=$((bad+1)); printf 'BADSIZE  %-20s %s\n' "${t}" "${dest}"; else ok=$((ok+1)); fi
done < "${MANIFEST}"
log "ok=${ok} bad_size=${bad} missing=${miss}"
