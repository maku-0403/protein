import os
import csv
import glob
from pathlib import Path
import pandas as pd

# ルートディレクトリ（あなたの環境に合わせて変更）
root_dir = input('Root directory: ')

# 全ての.csvファイルのフルパスを再帰的に取得
csv_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

process_count = 0
count = 0

#１つのCSVファイルごとに再帰的に繰り返し
for csv_path in csv_files:
    process_count += 1
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(csv_files)}")
    with open(csv_path) as f:
        reader = csv.reader(f)
        l = [row for row in reader]
        if len(l) != 0:
            if (len(l[0]) == 2) and (len(l) > 6):
                #search here
                if 'ELEC' in l[3][1]:
                    count += 1
print(count)