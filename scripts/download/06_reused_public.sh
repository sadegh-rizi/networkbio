#!/usr/bin/env bash
# 06_reused_public.sh — external datasets that Park et al. re-analysed (optional).
# Default is a dry run that lists the GEO supplementary files; pass --get to download.
#   GSE180759  Absinta 2021 snRNA-seq (MS chronic active lesions)
#   GSE279183  Lerma-Martin 2024 (subcortical MS lesions; Visium + snRNA; SuperSeries)
#   GSE208747  Alsema 2024 (WM lesion Visium)
#   GSE174647  Kaufmann 2022 (PMS grey-matter Visium)
# Schirmer 2019 (PRJNA544731) has raw reads only in SRA (~590 GB); its processed
# matrix is on the UCSC Cell Browser (dataset "ms"), which is not scripted here.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
get=0; [[ "${1:-}" == "--get" ]] && get=1
for g in GSE180759 GSE279183 GSE208747 GSE174647; do
  url="https://ftp.ncbi.nlm.nih.gov/geo/series/${g:0:${#g}-3}nnn/${g}/suppl/"
  log "${g}: ${url}"
  curl -fsL "${url}" | grep -o 'href="[^"?/][^"]*"' | sed 's/href="//;s/"$//' | grep -v '^http' | sed 's/^/    /'
  if (( get )); then
    for f in $(curl -fsL "${url}" | grep -o 'href="[^"?/][^"]*"' | sed 's/href="//;s/"$//' | grep -v '^http'); do
      fetch "${url}${f}" "${DATA_ROOT}/reused/${g}/${f}"
    done
  fi
done
