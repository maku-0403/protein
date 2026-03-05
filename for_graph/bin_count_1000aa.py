"""
1000AA未満の鎖を対象に、アミノ酸数（100AA刻み）× w≥10割合（%刻み）で
ビン分けし、テンプレート形式のCSVに出力するツール

X軸ビン（アミノ酸数）: 0-100AA, 100-200AA, ..., 900-1000AA  （1000未満のみ）
Y軸ビン（w≥10の割合）: 0-5%, 5-10%, 10-15%, 15-20%, 20-40%, 40-60%, 60-80%, 80-100%
"""

import os
import glob
import sys
import csv
import pandas as pd
from tqdm import tqdm


# ──────────────────────────────────────────────
# 1. ビン定義
# ──────────────────────────────────────────────

AA_BIN_EDGES  = list(range(0, 1001, 100))  # [0, 100, 200, ..., 1000]
AA_BIN_LABELS = [
    '0-100AA', '100-200AA', '200-300AA', '300-400AA', '400-500AA',
    '500-600AA', '600-700AA', '700-800AA', '800-900AA', '900-1000AA'
]

RATIO_BIN_EDGES  = [0, 5, 10, 15, 20, 40, 60, 80, 100]
RATIO_BIN_LABELS = [
    '0-5%', '5-10%', '10-15%', '15-20%',
    '20-40%', '40-60%', '60-80%', '80-100%'
]


def get_aa_bin(total_aa):
    """アミノ酸数のビンラベルを返す（1000以上はNone）"""
    if total_aa >= 1000:
        return None
    for i in range(len(AA_BIN_EDGES) - 1):
        if AA_BIN_EDGES[i] <= total_aa < AA_BIN_EDGES[i + 1]:
            return AA_BIN_LABELS[i]
    return None


def get_ratio_bin(ratio_pct):
    """割合(%)のビンラベルを返す（最終ビンは上限を含む）"""
    for i in range(len(RATIO_BIN_EDGES) - 1):
        lo = RATIO_BIN_EDGES[i]
        hi = RATIO_BIN_EDGES[i + 1]
        if i < len(RATIO_BIN_EDGES) - 2:
            if lo <= ratio_pct < hi:
                return RATIO_BIN_LABELS[i]
        else:  # 80-100%: 上限100%を含む
            if lo <= ratio_pct <= hi:
                return RATIO_BIN_LABELS[i]
    return None


# ──────────────────────────────────────────────
# 2. ユーザー入力
# ──────────────────────────────────────────────

print("=" * 50)
print(" 1000AA以下 ビン集計ツール")
print("=" * 50)

csv_dir = input("\nCSVフォルダのパスを入力: ").strip()
if not os.path.isdir(csv_dir):
    print(f"エラー: フォルダが見つかりません -> {csv_dir}")
    sys.exit(1)

out_path = input("出力CSVのパスを入力（例: /path/to/output.csv）: ").strip()
out_dir = os.path.dirname(out_path)
if out_dir and not os.path.exists(out_dir):
    os.makedirs(out_dir)

# フィルタリングCSV（任意）
filter_map = None
use_filter = input("\nフィルタリングCSV（pdb_chain_sp_primary_filtered.csv）を使用しますか？ [y/N]: ").strip().lower()
if use_filter == 'y':
    filter_csv_path = input("フィルタリングCSVのパスを入力: ").strip()
    if not os.path.isfile(filter_csv_path):
        print(f"エラー: ファイルが見つかりません -> {filter_csv_path}")
        sys.exit(1)
    filter_map = {}
    with open(filter_csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb = row['PDB'].upper()
            chain = row['CHAIN']
            filter_map.setdefault(pdb, set()).add(chain)
    print(f"フィルタリング対象: {len(filter_map)} PDB ID, "
          f"{sum(len(v) for v in filter_map.values())} CHAIN")


# ──────────────────────────────────────────────
# 3. CSVデータ読み込み・集計
# ──────────────────────────────────────────────

csv_files = glob.glob(os.path.join(csv_dir, "*.csv"))
if not csv_files:
    print(f"エラー: CSVファイルが見つかりません -> {csv_dir}")
    sys.exit(1)

if filter_map is not None:
    csv_files = [f for f in csv_files
                 if os.path.splitext(os.path.basename(f))[0].upper() in filter_map]
    print(f"\nフィルタリング後: {len(csv_files)} 件のCSVを処理中...")
else:
    print(f"\n{len(csv_files)} 件のCSVを処理中...")

# 集計テーブル初期化: {ratio_label: {aa_label: count}}
table = {r: {a: 0 for a in AA_BIN_LABELS} for r in RATIO_BIN_LABELS}

counted        = 0
skipped        = 0
skipped_ge1000 = 0

for fpath in tqdm(csv_files, unit="file"):
    try:
        with open(fpath, 'r', encoding='utf-8', errors='replace') as f:
            first_line = f.readline().strip()
        if first_line == 'NO_DATA':
            skipped += 1
            continue

        df = pd.read_csv(fpath, skiprows=8, header=0)
        df.columns = [c.strip() for c in df.columns]

        if 'amino_number' not in df.columns or 'w' not in df.columns:
            skipped += 1
            continue

        df['w'] = pd.to_numeric(df['w'], errors='coerce')
        df = df.dropna(subset=['w'])

        if len(df) == 0:
            skipped += 1
            continue

        fname = os.path.basename(fpath)

        if 'unit' in df.columns:
            df = df.dropna(subset=['unit'])
            if len(df) == 0:
                skipped += 1
                continue
            group_id = (df['unit'] != df['unit'].shift()).cumsum()
            groups = df.groupby(group_id, sort=False)
        else:
            groups = [(None, df)]

        for _, grp in groups:
            grp = grp.dropna(subset=['w'])
            if len(grp) == 0:
                continue

            total_aa = len(grp)
            if total_aa > 2500:
                continue

            unit_name = str(grp['unit'].iloc[0]) if 'unit' in grp.columns else 'N/A'

            if filter_map is not None:
                pdb_id = os.path.splitext(fname)[0].upper()
                if unit_name not in filter_map.get(pdb_id, set()):
                    continue

            # 1000AA以上はこの集計対象外
            if total_aa >= 1000:
                skipped_ge1000 += 1
                continue

            ratio     = (grp['w'] >= 10).sum() / len(grp)
            ratio_pct = ratio * 100

            aa_label    = get_aa_bin(total_aa)
            ratio_label = get_ratio_bin(ratio_pct)

            if aa_label is not None and ratio_label is not None:
                table[ratio_label][aa_label] += 1
                counted += 1

    except Exception:
        skipped += 1

print(f"\n集計完了: {counted} 件カウント"
      f" / 1000AA以上スキップ: {skipped_ge1000} 件"
      f" / ファイルスキップ: {skipped} 件")

if counted == 0:
    print("警告: 集計対象データが0件でした。")


# ──────────────────────────────────────────────
# 4. CSV出力
# ──────────────────────────────────────────────

with open(out_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow([''] + AA_BIN_LABELS)
    for r_label in RATIO_BIN_LABELS:
        row = [r_label] + [table[r_label][a] for a in AA_BIN_LABELS]
        writer.writerow(row)

print(f"CSV出力完了: {out_path}")
