#!/usr/bin/env bash
# 02_atac_fragments.sh — snATAC fragments.tsv.gz for GSE297690 (4 libraries, ~6.2 GB).
# Needed for Signac/ArchR re-analysis; the peak .h5 files from 01 are enough for a
# first look.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_space 8
fetch_tier geo_atac_fragments
