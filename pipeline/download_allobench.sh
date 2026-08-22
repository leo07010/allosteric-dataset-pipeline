#!/bin/bash
# 下載 AlloBench 所需的 PDB 結構。RCSB 有速率限制,用 8 個並行 + 每個間隔,
# 並跳過已存在的檔(可中斷後重跑)。
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$REPO_ROOT"
OUT="${PDB_DIR:-data/allobench_pdbs}"
mkdir -p "$OUT"
python3 -c "
import json
import os
d=json.load(open(os.environ.get('FILTERED_JSON','metadata/allobench_filtered.json')))
print('\n'.join(d['pdbs']))" > /tmp/pdb_list.txt
TOTAL=$(wc -l < /tmp/pdb_list.txt)
echo "需處理 $TOTAL 個 PDB"
cat /tmp/pdb_list.txt | xargs -P 8 -I{} sh -c '
  f='"$OUT"'/{}.pdb
  [ -s "$f" ] && exit 0
  # 先找本地既有的
  curl -sS --max-time 45 -o "$f" "https://files.rcsb.org/download/{}.pdb" || rm -f "$f"
  [ -s "$f" ] || rm -f "$f"
  sleep 0.15
'
echo "完成: $(ls $OUT/*.pdb 2>/dev/null | wc -l) / $TOTAL"
