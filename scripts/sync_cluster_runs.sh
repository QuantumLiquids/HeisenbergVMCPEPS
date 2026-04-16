#!/bin/bash
# Sync cluster run directories to local, excluding large binary data.
# Safe: never deletes local-only files (e.g. run_history.json).
#
# Usage: ./scripts/sync_cluster_runs.sh [subdirectory]
#   No args: sync all run/*
#   With arg: sync only run/<subdirectory>  (e.g. "8x8J2=0.5D8")

set -euo pipefail

REMOTE="wanghx@10.20.26.130"
PROXY_CMD="nc -X 5 -x 127.0.0.1:1080 %h %p"
REMOTE_BASE="/share/home/wanghx/HeisenbergVMCPEPS/run/"
LOCAL_BASE="/Users/wanghaoxin/GitHub/HeisenbergVMCPEPS/run/"

# What to exclude (large binary data we don't need locally)
EXCLUDES=(
    '*.qlten'           # tensor files (MB each, thousands per run)
    'configuration*'    # MC configuration snapshots
    'tpsfinal/'         # full TPS directory
    'tpslowest/'        # lowest-energy TPS snapshots
    'tps_checkpoint*'   # checkpoint TPS directories
    'peps/'             # raw PEPS tensors
    'measure/'          # measurement output tensors
    '*.core'            # core dump files
    '*.lz4'             # compressed core dumps
)

# Build exclude flags
EXCLUDE_FLAGS=""
for pat in "${EXCLUDES[@]}"; do
    EXCLUDE_FLAGS="$EXCLUDE_FLAGS --exclude=$pat"
done

# Determine what to sync
if [ $# -gt 0 ]; then
    SUBDIR="$1/"
    echo "=== Syncing run/$1 ==="
else
    SUBDIR=""
    echo "=== Syncing all run/ directories ==="
fi

# rsync: archive mode, compress, show progress
# NO --delete: preserves local-only files (run_history.json, notes, etc.)
# --update: skip files newer on local (protects local edits)
rsync -avz --progress \
    --update \
    -e "ssh -o ProxyCommand=\"$PROXY_CMD\"" \
    $EXCLUDE_FLAGS \
    "${REMOTE}:${REMOTE_BASE}${SUBDIR}" \
    "${LOCAL_BASE}${SUBDIR}"

echo "=== Sync complete ==="
echo "Note: local-only files (run_history.json etc.) are preserved."
echo "      Large files (*.qlten, configurations, tpsfinal/) are excluded."
