import csv
import glob
import os

# ルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/Users/kuniimahan/Desktop/school/研究関係/千葉研/protein-1/for_csv"

# 全ての.csvファイルのフルパスを再帰的に取得
csv_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

process_count = 0
all_list = list()

#１つのCSVファイルごとに再帰的に繰り返し
for csv_path in csv_files:
    process_count += 1
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(csv_files)}")
    with open(csv_path) as f:
        reader = csv.reader(f)
        l = [row for row in reader]

with open(csv_path,'w',newline="") as f:
        writer = csv.writer(f)
        for i in range(0,len(l)):
            if len(l[i]) != 0:
                writer.writerow(l[i])