#!/usr/bin/env bash
#
# Fetch the computational mesh needed to run the WINDSOR case:
#   - MESH/c1g1.cgns  (~1.2 GB)  ~6.3 M cells
#
# The mesh is hosted outside the git repository because of its size.
# Once downloaded, the case directories are self-contained and can be
# run with the standard `code_saturne run` workflow.
#
# Usage:
#   ./fetch_data.sh
#
# The download URL is configured via the WINDSOR_DATA_URL environment
# variable, or by editing the DATA_URL variable below.

set -euo pipefail

DATA_URL="${WINDSOR_DATA_URL:-<REPLACE_WITH_HOSTED_URL>}"

cd "$(dirname "$0")"

if [ "$DATA_URL" = "<REPLACE_WITH_HOSTED_URL>" ]; then
    echo "ERROR: download URL not set." >&2
    echo "  Either edit DATA_URL in this script or export WINDSOR_DATA_URL." >&2
    exit 1
fi

fetch() {
    local rel_path="$1"
    local sha256="${2:-}"

    if [ -f "$rel_path" ]; then
        echo "[skip] $rel_path already present"
        return
    fi

    mkdir -p "$(dirname "$rel_path")"
    echo "[get ] $rel_path"
    curl -fL --progress-bar "$DATA_URL/$rel_path" -o "$rel_path"

    if [ -n "$sha256" ]; then
        echo "[chk ] verifying sha256"
        echo "$sha256  $rel_path" | sha256sum -c -
    fi
}

fetch "MESH/c1g1.cgns"

echo "Done. You can now run:"
echo "  code_saturne run --case RANS_G1"
echo "  code_saturne run --case DDES_G1"
