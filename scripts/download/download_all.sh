#!/usr/bin/env bash
# download_all.sh — convenience wrapper. Default: processed files only (~4 GB).
#   bash download_all.sh              processed + metadata + code
#   bash download_all.sh --medium     + ATAC fragments + raw LC-MS (~15 GB more)
#   Raw sequencing reads are never fetched here; use 04_sra_fastq.sh explicitly.
set -euo pipefail
here="$(dirname "${BASH_SOURCE[0]}")"
bash "${here}/01_processed.sh"
if [[ "${1:-}" == "--medium" ]]; then
  bash "${here}/02_atac_fragments.sh"
  bash "${here}/03_mw_raw.sh"
fi
bash "${here}/05_massive_proteomics.sh" || echo "MassIVE step did not complete (dataset private?)"
bash "${here}/verify.sh"
