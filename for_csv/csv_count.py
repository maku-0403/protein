import csv
import os
import glob
from pathlib import Path
import pandas as pd

# sphereルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/Volumes/pdb_res/CIF/allcsv_to_retouch/0.3"

# 全ての.csvファイルのフルパスを再帰的に取得
csv_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

process_count = 0

#１つのCSVファイルごとに再帰的に繰り返し
for csv_path in csv_files:
    process_count += 1
print(process_count)