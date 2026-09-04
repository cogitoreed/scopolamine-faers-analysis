#!/bin/bash
## FAERS 并行下载主控（8 并发）
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$DIR/data"; mkdir -p "$DATA_DIR"
LOG="$DATA_DIR/parallel_download.log"; : > "$LOG"
QFILE="$DATA_DIR/quarters.txt"; : > "$QFILE"
for y in $(seq 2012 2026); do
  for q in 1 2 3 4; do
    [ "$y" -eq 2012 ] && [ "$q" -ne 4 ] && continue
    [ "$y" -eq 2026 ] && [ "$q" -gt 2 ] && continue
    printf '%dQ%d\n' "$y" "$q" >> "$QFILE"
  done
done
echo "共 $(wc -l < "$QFILE") 个季度，启动 8 并发..."
xargs -P 8 -I Q "$DIR/02a_dl_one.sh" Q < "$QFILE"
echo "=== 全部结束 成功/跳过: $(grep -c '^\[ok\]\|^\[skip\]' "$LOG") 失败: $(grep -c '^\[FAIL\]' "$LOG") ===" >> "$LOG"
