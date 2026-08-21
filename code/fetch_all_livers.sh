#!/usr/bin/env bash
# Sequential fetch of the remaining GSE310392 liver bundles.
# Sequential, not parallel: avoids hammering NCBI, and avoids a repeat of the
# quota race that silently truncated the first two downloads.
set -u
cd /workspace
while read -r gsm fn; do
  [ -z "$gsm" ] && continue
  echo "=== $(date -u +%FT%TZ) $gsm $fn ==="
  bash code/fetch_xenium_bundle.sh "$gsm" "$fn" /workspace/data/raw \
    && echo "  OK" || echo "  FAILED $gsm"
  du -sh /workspace | awk '{print "  volume now: "$1}'
done < "${1:-logs/sample_list.txt}"
echo "=== ALL DONE $(date -u +%FT%TZ) ==="
ls -d /workspace/data/raw/*/ | sed 's|.*/raw/||'
