#!/bin/bash
## 下载单个季度: 02a_dl_one.sh 2013Q1
set -u
BASE="https://fis.fda.gov/content/Exports/faers_ascii"
DATA_DIR="$(cd "$(dirname "$0")" && pwd)/data"
LOG="$DATA_DIR/parallel_download.log"
tag="$1"
tagl=$(echo "$tag" | tr 'A-Z' 'a-z')
dest="$DATA_DIR/faers_ascii_${tagl}.zip"
sz=$(stat -f%z "$dest" 2>/dev/null || echo 0)
if [ "$sz" -gt 1000000 ]; then echo "[skip] $tag ($((sz/1048576))MB)" >> "$LOG"; exit 0; fi
for t in "$tag" "$tagl"; do
  if curl -sfL -C - --retry 4 --retry-delay 3 --max-time 14400 --speed-limit 2000 --speed-time 180 \
       -o "$dest" "${BASE}_${t}.zip" 2>>"$DATA_DIR/curl_err.log"; then
    echo "[ok]   $tag ($(( $(stat -f%z "$dest") /1048576 ))MB)" >> "$LOG"; exit 0
  fi
done
echo "[FAIL] $tag" >> "$LOG"; rm -f "$dest"; exit 1
