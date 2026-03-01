"""
CSV 分子量可視化ツール
- 散布図（密度なし）
- 散布図（密度あり・KDE色分け）
- 3D棒グラフ（X×Yビンごとの個数をZ軸）
"""

import os
import glob
import sys
import csv
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm

matplotlib.use('Agg')


# ──────────────────────────────────────────────
# 1. ユーザー入力
# ──────────────────────────────────────────────

print("=" * 50)
print(" CSV 分子量可視化ツール")
print("=" * 50)

print("\n【プロット種別を選択してください】")
print("  1: 散布図（密度なし）")
print("  2: 散布図（密度あり・KDE色分け）")
print("  3: 3D棒グラフ（個数軸）")

while True:
    choice = input("\n番号を入力 [1/2/3]: ").strip()
    if choice in ("1", "2", "3"):
        break
    print("  ※ 1, 2, 3 のいずれかを入力してください。")

csv_dir = input("\nCSVフォルダのパスを入力: ").strip()
if not os.path.isdir(csv_dir):
    print(f"エラー: フォルダが見つかりません -> {csv_dir}")
    sys.exit(1)

out_path = input("出力ファイルのパスを入力（例: /path/to/output.png）: ").strip()
out_dir = os.path.dirname(out_path)
if out_dir and not os.path.exists(out_dir):
    os.makedirs(out_dir)

log_path = os.path.splitext(out_path)[0] + ".log.csv"

# フィルタリングCSV（任意）
# pdb_chain_sp_primary_filtered.csv を指定すると、
# 対応する {PDBid}.csv の指定 CHAIN のみをプロット対象にする
filter_map = None  # {PDB_ID(大文字): set(CHAIN)}
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
# 2. CSVデータ読み込み
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

x_values    = []
y_values    = []
log_records = []  # (filename, unit, amino_count, ratio)
skipped     = 0

for fpath in tqdm(csv_files, unit="file"):
    try:
        # NO_DATA ファイルをスキップ
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

        # unit 列があれば行順で切り替わりを検出してグループ化
        # NaN unit 行はユニット不明のため除外する
        if 'unit' in df.columns:
            df = df.dropna(subset=['unit'])          # ← NaN unit 行を除外
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

            # CHAINフィルタリング
            if filter_map is not None:
                pdb_id = os.path.splitext(fname)[0].upper()
                if unit_name not in filter_map.get(pdb_id, set()):
                    continue

            w_vals = grp['w']
            ratio = (w_vals >= 10).sum() / len(w_vals)
            x_values.append(total_aa)
            y_values.append(ratio)
            log_records.append((fname, unit_name, total_aa, ratio))

    except Exception:
        skipped += 1

processed = len(x_values)
print(f"処理済み: {processed} プロット点 / スキップファイル: {skipped} 件")

if processed == 0:
    print("エラー: プロット可能なデータがありません。")
    sys.exit(1)

# ──────────────────────────────────────────────
# 3. ログCSV出力
# ──────────────────────────────────────────────

with open(log_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['filename', 'unit', 'amino_count', 'ratio_w_ge10'])
    writer.writerows(log_records)
print(f"ログCSVを保存しました: {log_path}")

# ──────────────────────────────────────────────
# 4. プロット
# ──────────────────────────────────────────────

x_arr = np.array(x_values)
y_arr = np.array(y_values)

TITLE_BASE = f"amino acid count vs proportion of w ≥ 10  (N={processed})"

# ---------- 選択肢 1: 散布図（密度なし）----------
if choice == "1":
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(x_arr, y_arr, alpha=0.4, s=10, color='steelblue', edgecolors='none')
    ax.set_xlabel("Total amino acid count", fontsize=13)
    ax.set_ylabel("Proportion of w ≥ 10", fontsize=13)
    ax.set_title("Scatter plot: " + TITLE_BASE, fontsize=13)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle='--', alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)

# ---------- 選択肢 2: 散布図（密度あり）----------
elif choice == "2":
    from scipy.stats import gaussian_kde

    xy = np.vstack([x_arr, y_arr])
    density = gaussian_kde(xy)(xy)

    sort_idx = np.argsort(density)
    x_s, y_s, d_s = x_arr[sort_idx], y_arr[sort_idx], density[sort_idx]

    fig, ax = plt.subplots(figsize=(10, 7))
    sc = ax.scatter(x_s, y_s, c=d_s, cmap='jet', alpha=0.6, s=10, edgecolors='none')
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Density (KDE)", fontsize=11)
    ax.set_xlabel("Total amino acid count", fontsize=13)
    ax.set_ylabel("Proportion of w ≥ 10", fontsize=13)
    ax.set_title("Scatter plot (KDE density): " + TITLE_BASE, fontsize=13)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle='--', alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)

# ---------- 選択肢 3: 3D棒グラフ ----------
elif choice == "3":
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    x_bins, y_bins = 60, 40
    hist, x_edges, y_edges = np.histogram2d(x_arr, y_arr, bins=[x_bins, y_bins])

    xpos = (x_edges[:-1] + x_edges[1:]) / 2
    ypos = (y_edges[:-1] + y_edges[1:]) / 2
    dx = (x_edges[1] - x_edges[0]) * 0.8
    dy = (y_edges[1] - y_edges[0]) * 0.8

    xg, yg = np.meshgrid(xpos, ypos, indexing='ij')
    xf, yf = xg.ravel(), yg.ravel()
    zf = np.zeros_like(xf)
    dz = hist.ravel()

    mask = dz > 0
    xf, yf, zf, dz = xf[mask], yf[mask], zf[mask], dz[mask]

    norm = plt.Normalize(vmin=dz.min(), vmax=dz.max())
    colors = plt.cm.jet(norm(dz))

    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(xf, yf, zf, dx, dy, dz, color=colors, zsort='average', alpha=0.85)

    mappable = plt.cm.ScalarMappable(cmap='jet', norm=norm)
    mappable.set_array(dz)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.5, pad=0.1)
    cbar.set_label("Count", fontsize=11)

    ax.set_xlabel("Total amino acid count", fontsize=11, labelpad=10)
    ax.set_ylabel("Proportion of w ≥ 10", fontsize=11, labelpad=10)
    ax.set_zlabel("Count", fontsize=11, labelpad=8)
    ax.set_title("3D histogram: " + TITLE_BASE, fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)

print(f"保存しました: {out_path}")
