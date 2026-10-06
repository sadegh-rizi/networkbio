#!/usr/bin/env bash
# 04_sra_fastq.sh — raw reads for the Park GEO series.
#
# Usage: bash 04_sra_fastq.sh [--method ena|sra] [--dry-run] SERIES [SERIES ...]
#   SERIES: GSE297192 (bulk RNA, ~50 GB)  GSE297365 (scRNA, ~38 GB)
#           GSE297690 (snATAC, ~38 GB)    GSE251839 (WGBS, >880 GB)
#   --method ena (default): wget the gzipped FASTQ from ENA. No extra software.
#   --method sra: SRA Toolkit prefetch + fasterq-dump (needs sra-tools on PATH).
#            Use this for 10x data if you want the index reads
#            (--include-technical) for Cell Ranger, and for runs ENA has not mirrored.
#
# Run list and sizes are in sra_runs.tsv (from NCBI SRA runinfo + ENA, 2026-09-29).
# The 6 WGBS runs SRR33556538-43 reported size 0 in SRA and no ENA FASTQ on that
# date; they are attempted with the sra method and skipped with a warning otherwise.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

method=ena; dry=0; series=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --method) method="$2"; shift 2;;
    --dry-run) dry=1; shift;;
    GSE*) series+=("$1"); shift;;
    *) echo "unknown argument: $1"; exit 2;;
  esac
done
[[ ${#series[@]} -eq 0 ]] && { sed -n '2,16p' "$0"; exit 2; }

declare -A GB=([GSE297192]=55 [GSE297365]=45 [GSE297690]=45 [GSE251839]=1000)
total=0; for s in "${series[@]}"; do total=$(( total + ${GB[$s]:-50} )); done
(( dry )) || need_space "${total}"

for s in "${series[@]}"; do
  out="${DATA_ROOT}/Park/SRA/${s}"; mkdir -p "${out}"
  # columns: series run gsm sample_title condition strategy spots sra_size_MB ena_fastq_ftp ena_fastq_bytes
  while IFS=$'\t' read -r -u 3 ser run gsm title cond strat _spots mb ftp bytes; do
    [[ "${ser}" != "${s}" ]] && continue
    [[ "${ftp}" == "-" ]] && ftp=""; [[ "${bytes}" == "-" ]] && bytes=""
    log "${s} ${run} ${gsm} ${title} (${cond}; SRA ${mb} MB)"
    if (( dry )); then continue; fi
    if [[ "${method}" == "ena" && -n "${ftp}" ]]; then
      IFS=';' read -ra urls <<< "${ftp}"; IFS=';' read -ra sizes <<< "${bytes}"
      for i in "${!urls[@]}"; do
        f="${urls[$i]##*/}"
        fetch "https://${urls[$i]}" "${out}/${gsm}_${title}/${f}" "${sizes[$i]:-}"
      done
    else
      if ! have prefetch || ! have fasterq-dump; then
        log "WARN: ${run}: no ENA FASTQ and sra-tools not installed; skipped"; continue
      fi
      if [[ "${mb}" == "0" ]]; then log "WARN: ${run} reports 0 MB in SRA; trying anyway"; fi
      d="${out}/${gsm}_${title}"; mkdir -p "${d}"
      prefetch --max-size u --output-directory "${d}" "${run}" || { log "WARN: prefetch failed for ${run}"; continue; }
      extra=(); [[ "${strat}" != "Bisulfite-Seq" && "${s}" != "GSE297192" ]] && extra=(--include-technical)
      fasterq-dump --split-files "${extra[@]}" --threads "${THREADS:-4}" --outdir "${d}" "${d}/${run}/${run}.sra"
      gzip -f "${d}/${run}"*.fastq
      rm -rf "${d:?}/${run}"
    fi
  done 3< "${RUNS_TSV}"
done
log "done"
