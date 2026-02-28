"""
CSV 分子量可視化ツール
- 散布図（密度なし）
- 散布図（密度あり・KDE色分け）
- 3D棒グラフ（X×Yビンごとの個数をZ軸）
- インタラクティブ散布図（plotly・HTML出力、ジッターあり）
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
print("  4: インタラクティブ散布図（HTML・ホバー情報付き）")

while True:
    choice = input("\n番号を入力 [1/2/3/4]: ").strip()
    if choice in ("1", "2", "3", "4"):
        break
    print("  ※ 1, 2, 3, 4 のいずれかを入力してください。")

csv_dir = input("\nCSVフォルダのパスを入力: ").strip()
if not os.path.isdir(csv_dir):
    print(f"エラー: フォルダが見つかりません -> {csv_dir}")
    sys.exit(1)

default_ext = ".html" if choice == "4" else ".png"
out_path = input(f"出力ファイルのパスを入力（例: /path/to/output{default_ext}）: ").strip()
out_dir = os.path.dirname(out_path)
if out_dir and not os.path.exists(out_dir):
    os.makedirs(out_dir)

log_path = os.path.splitext(out_path)[0] + ".log.csv"

# ──────────────────────────────────────────────
# 2. CSVデータ読み込み
# ──────────────────────────────────────────────

csv_files = glob.glob(os.path.join(csv_dir, "*.csv"))
if not csv_files:
    print(f"エラー: CSVファイルが見つかりません -> {csv_dir}")
    sys.exit(1)

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

# ---------- 選択肢 4: インタラクティブ散布図（plotly）----------
elif choice == "4":
    import plotly.graph_objects as go

    filenames  = [r[0] for r in log_records]
    unit_names = [r[1] for r in log_records]
    aa_counts  = [r[2] for r in log_records]
    ratios     = [r[3] for r in log_records]

    # X軸（アミノ酸数）にジッターを加えて重なりを緩和
    jitter_scale = max(aa_counts) * 0.003
    rng = np.random.default_rng(seed=42)
    x_jitter = x_arr + rng.uniform(-jitter_scale, jitter_scale, size=len(x_arr))

    hover_text = [
        f"<b>File:</b> {fn}<br>"
        f"<b>Unit:</b> {un}<br>"
        f"<b>Amino acids:</b> {aa}<br>"
        f"<b>Ratio (w≥10):</b> {rt:.4f}"
        for fn, un, aa, rt in zip(filenames, unit_names, aa_counts, ratios)
    ]

    fig = go.Figure(go.Scatter(
        x=x_jitter,
        y=y_arr,
        mode='markers',
        marker=dict(size=5, color='steelblue', opacity=0.5),
        text=hover_text,
        hovertemplate="%{text}<extra></extra>",
    ))

    fig.update_layout(
        title=dict(text="Interactive scatter: " + TITLE_BASE, font=dict(size=15)),
        xaxis=dict(title="Total amino acid count", showgrid=True, gridcolor='lightgrey'),
        yaxis=dict(title="Proportion of w ≥ 10", showgrid=True, gridcolor='lightgrey',
                   range=[-0.05, 1.05]),
        plot_bgcolor='white',
        hoverlabel=dict(bgcolor='white', font_size=13),
        width=1000,
        height=700,
    )

    fig.write_html(out_path, include_plotlyjs='cdn')

print(f"保存しました: {out_path}")
