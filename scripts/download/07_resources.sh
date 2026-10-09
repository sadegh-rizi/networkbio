#!/usr/bin/env bash
# Download the approved HGNC quarterly release; never replace existing resources.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEST="${DATA_ROOT:-${ROOT}/data}/resources/hgnc"
PYTHON="${ROOT}/.venv/bin/python"
[[ -x "${PYTHON}" ]] || { printf 'Missing project Python: %s\n' "${PYTHON}" >&2; exit 1; }
command -v curl >/dev/null || { printf 'curl is required.\n' >&2; exit 1; }
mkdir -p "${DEST}"
BASE="https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv"
selected=""
downloaded=false

# A cached fallback remains the selected resource on subsequent runs.
for release in 2026-10-06 2026-07-07; do
    candidate="${DEST}/hgnc_complete_set_${release}.txt"
    if [[ -e "${candidate}" ]]; then
        selected="${candidate}"
        break
    fi
done
if [[ -z "${selected}" ]]; then
    temp="$(mktemp "${DEST}/.hgnc.XXXXXX")"
    trap 'rm -f "${temp}"' EXIT
    for release in 2026-10-06 2026-07-07; do
        candidate="${DEST}/hgnc_complete_set_${release}.txt"
        url="${BASE}/hgnc_complete_set_${release}.txt"
        status="$(curl -sSL --retry 3 --connect-timeout 30 --max-time 300 -w '%{http_code}' -o "${temp}" "${url}")"
        case "${status}" in
            200)
                # Check schema before publishing the new file (a 200 HTML page is not HGNC).
                "${PYTHON}" - "${temp}" <<'PY'
import sys
import pandas as pd
table = pd.read_csv(sys.argv[1], sep="\t", nrows=5)
if table.empty or not {"ensembl_gene_id", "symbol", "status"}.issubset(table.columns):
    raise ValueError("Downloaded file is not the expected HGNC TSV")
PY
                ln "${temp}" "${candidate}"  # fails rather than overwriting an existing resource
                selected="${candidate}"
                downloaded=true
                break
                ;;
            404|410) printf 'HGNC release %s unavailable (HTTP %s).\n' "${release}" "${status}" ;;
            *) printf 'HGNC download failed (HTTP %s): %s\n' "${status}" "${url}" >&2; exit 1 ;;
        esac
    done
fi
[[ -n "${selected}" ]] || { printf 'Neither approved HGNC release is available.\n' >&2; exit 1; }
"${PYTHON}" - "${selected}" "${BASE}" "${downloaded}" <<'PY'
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

path = Path(sys.argv[1])
sidecar = path.with_suffix(".provenance.json")
digest = hashlib.sha256(path.read_bytes()).hexdigest()
if sys.argv[3] == "true":
    table = pd.read_csv(path, sep="\t", usecols=["ensembl_gene_id", "symbol", "status"])
    metadata = {"resource": "HGNC complete set", "release": path.stem.removeprefix("hgnc_complete_set_"),
                "url": f"{sys.argv[2]}/{path.name}", "downloaded_utc": datetime.now(timezone.utc).isoformat(),
                "sha256": digest, "row_count": len(table)}
    with sidecar.open("x") as handle:
        json.dump(metadata, handle, indent=2)
elif not sidecar.exists():
    raise ValueError(f"Existing HGNC resource has no download provenance: {path}; inspect it before proceeding")
else:
    metadata = json.loads(sidecar.read_text())
    if metadata["sha256"] != digest:
        raise ValueError(f"HGNC checksum mismatch: {path}")
print(f"HGNC resource ready: {path} (sha256 {digest})")
PY
