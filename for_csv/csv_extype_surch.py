import os
import csv
import glob
from pathlib import Path

# ルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/Volumes/pdb_res/CIF/cif_to_csv/all_csv_cosw"

# 全ての.csvファイルのフルパスを再帰的に取得
csv_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

process_count = 0
count = 0

exp_list = list()

#１つのCSVファイルごとに再帰的に繰り返し
for csv_path in csv_files:
    process_count += 1
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(csv_files)}")
    with open(csv_path) as f:
        reader = csv.reader(f)
        l = [row for row in reader]
        if len(l) != 0:
            if (len(l[0]) == 2) and (len(l) > 2):
                #search here
                if len(exp_list) == 0:
                    exp_list.append([l[3][1],1])
                else:
                    for i in range(0,len(exp_list)):
                        if exp_list[i][0] == l[3][1]:
                            exp_list[i][1] += 1
                            break
                        if i == len(exp_list)-1:
                            exp_list.append([l[3][1],1])

with open('exp_list.csv', 'w') as f:
    writer = csv.writer(f)
    for i in range(0,len(exp_list)):
        writer.writerow([exp_list[i][0],exp_list[i][1]])
print(exp_list)