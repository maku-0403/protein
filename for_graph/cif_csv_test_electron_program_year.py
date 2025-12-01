import csv
import glob
import os
from collections import Counter, defaultdict

# mmCIFのCSVが置かれているルートディレクトリ（必要に応じて変更）
root_dir = input("CSV dir: ")

# 出力先（必要に応じて変更）
output_dir = input("Out dir: ")
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "program_count_by_year.csv")

# 全ての.csvファイルを再帰的に取得
csv_files = glob.glob(os.path.join(root_dir, "**", "*.csv"), recursive=True)

# year -> Counter(program_name -> count)
counts_by_year: dict[str, Counter] = defaultdict(Counter)
program_names = set()

for idx, csv_path in enumerate(csv_files, start=1):
    if idx % 100 == 0:
        print(f"処理中: {idx} / {len(csv_files)}")

    with open(csv_path) as f:
        reader = csv.reader(f)
        rows = [row for row in reader]

    # 想定外のフォーマットはスキップ
    if not rows or len(rows[0]) != 2 or len(rows) < 5:
        continue

    # 実験手法フィルタ（X-RAY/NEUTRONのみ）
    experiment = rows[3][1] if len(rows[3]) > 1 else ""
    if "ELECTRON" not in experiment:
        continue

    # 年取得
    date_str = rows[2][1] if len(rows[2]) > 1 else ""
    if len(date_str) < 4:
        continue
    year = date_str[:4]

    # プログラム名取得（row 4 の2列目以降）
    program_row = rows[4]
    programs = {
        p.strip()
        for p in program_row[1:]
    }
    if not programs:
        continue

    for program in programs:
        counts_by_year[year][program] += 1
        program_names.add(program)

# 並び順を固定
sorted_years = sorted(counts_by_year.keys())
sorted_programs = sorted(program_names)

with open(output_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["year"] + sorted_programs)
    for year in sorted_years:
        row = [counts_by_year[year].get(program, 0) for program in sorted_programs]
        writer.writerow([year] + row)

print(f"Saved: {output_path}")
