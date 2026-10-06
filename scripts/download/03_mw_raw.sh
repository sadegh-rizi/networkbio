#!/usr/bin/env bash
# 03_mw_raw.sh — raw LC-MS files (Thermo .raw in zip) for ST003328/30/31/32 (~8.7 GB).
# Only needed to re-process spectra; 01_processed.sh already gets the peak tables.
# ST003328 lists a second zip (ST003328_Rawfiles.zip) with the same byte size as
# ST003328_DS1-93_Rawdata.zip; it is skipped unless you pass --with-duplicate.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_space 10
fetch_tier mw_raw
if [[ "${1:-}" == "--with-duplicate" ]]; then need_space 5; fetch_tier mw_raw_duplicate; fi
