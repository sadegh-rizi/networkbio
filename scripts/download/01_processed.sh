#!/usr/bin/env bash
# 01_processed.sh — small, analysis-ready files (about 4 GB in total).
#   * GEO SOFT + series-matrix metadata for the 4 Park series
#   * Park processed files, per sample (bulk count matrices, scRNA .h5,
#     snATAC peak .h5 + fragment index, WGBS CpG tables / Bismark .cov)
#   * Ionescu Metabolomics Workbench processed data (REST JSON, mwTab, factors)
#   * Code archives from Zenodo
# ATAC fragment files (6.2 GB) are in 02_atac_fragments.sh.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_space 6
for tier in meta geo_processed mw_processed code; do fetch_tier "${tier}"; done
log "done: processed tier -> ${DATA_ROOT}"
