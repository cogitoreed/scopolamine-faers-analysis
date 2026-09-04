#!/bin/bash
###############################################################################
## FAERS / AEMS 季度数据批量下载器
##
## 用法：在终端执行
##   cd "/Users/xiaocaixu/WorkBuddy/学位自救计划/FAERS分析"
##   chmod +x 01_下载数据.sh
##   ./01_下载数据.sh
##
## 说明：
##   * FDA 官方目前只提供 2012Q4 至 2026Q2 共 55 个季度，更早的季度返回 500
##   * 服务器文件名大小写混乱（有的 2025Q4，有的 2025q3），脚本会自动两种都试
##   * 支持断点续传，中断后重新执行即可，已完成的会跳过
##   * 全量约 3.5 GB，视网速约需 20-60 分钟，建议挂在后台跑
##
## 后台运行（关掉终端也不会断）：
##   nohup ./01_下载数据.sh > download.log 2>&1 &
##   查看进度：tail -f download.log
###############################################################################

set -u
DATA_DIR="$(cd "$(dirname "$0")" && pwd)/data"
mkdir -p "$DATA_DIR"
BASE="https://fis.fda.gov/content/Exports/faers_ascii"

ok=0; fail=0; skip=0; failed_list=()

for year in $(seq 2012 2026); do
  for q in 1 2 3 4; do
    # 2012 只有 Q4；2026 目前只有 Q1、Q2
    if [ "$year" -eq 2012 ] && [ "$q" -ne 4 ]; then continue; fi
    if [ "$year" -eq 2026 ] && [ "$q" -gt 2 ]; then continue; fi

    tag="$(printf '%dQ%d' "$year" "$q")"
    tagl="$(printf '%dq%d' "$year" "$q")"
    dest="$DATA_DIR/faers_ascii_${tagl}.zip"

    # 已完成且大小合理则跳过
    if [ -f "$dest" ]; then
      sz=$(stat -f%z "$dest" 2>/dev/null || echo 0)
      if [ "$sz" -gt 1000000 ]; then
        echo "[跳过] $tag 已存在 ($(( sz / 1048576 )) MB)"
        skip=$((skip+1)); continue
      fi
    fi

    # 大小写两种都试，带断点续传
    got=0
    for t in "$tag" "$tagl"; do
      code=$(curl -sIL -o /dev/null -w '%{http_code}' --max-time 20 \
             "${BASE}_${t}.zip")
      if [ "$code" = "200" ]; then
        echo "[下载] $t ..."
        if curl -fL -C - --retry 3 --retry-delay 3 --max-time 1800 \
                -o "$dest" "${BASE}_${t}.zip"; then
          got=1; break
        fi
      fi
    done

    if [ "$got" -eq 1 ]; then
      sz=$(stat -f%z "$dest" 2>/dev/null || echo 0)
      echo "       完成 $(( sz / 1048576 )) MB"
      ok=$((ok+1))
    else
      echo "[失败] $tag —— 两种大小写均不可用"
      failed_list+=("$tag"); fail=$((fail+1))
      rm -f "$dest"
    fi
  done
done

echo ""
echo "==================== 下载汇总 ===================="
echo "成功: $ok    跳过: $skip    失败: $fail"
echo "文件目录: $DATA_DIR"
if [ "$fail" -gt 0 ]; then
  echo "失败季度: ${failed_list[*]}"
  echo "→ 可打开 https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html 手动补下"
fi
echo "=================================================="
