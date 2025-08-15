import matplotlib.pyplot as plt
import csv
import os
import glob

# mmCIFルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/Volumes/pdb_res/CIF/cif_to_csv/all_csv_temperature/sphere"
res_dir = "/Volumes/pdb_res/CIF/cif_to_csv/all_csv_cosw/"

# 全ての.cifファイルのフルパスを再帰的に取得
csv_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

x = list()
y = list()

process_count = 0

for csv_path in csv_files:
    process_count += 1
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(csv_files)}")
    if process_count == 1000:
        break
    with open(csv_path) as f:
            reader = csv.reader(f)
            l = [row for row in reader]
    if len(l) != 0:
        if (len(l[0]) == 2) and (len(l) > 2):
            try:
                temp = int(l[2][1])
                sum = 0
                for i in range(2,len(l)):
                    sum += float(l[i][3])
                avg = sum / (len(l)-2)
                if avg != 0:
                    x.append(avg)
                    file_name = os.path.splitext(os.path.basename(path))[0]
                    y_path = res_dir+file_name+'.csv'
                    with open(y_path) as f:
                        reader = csv.reader(f)
                        l = [row for row in reader]
                        resolution = l[1][1]
                        y.append(resolution)
            except:
                continue

# 散布図の作成
plt.scatter(x, y)

# 軸ラベル（任意）
plt.xlabel("Resolution")
plt.ylabel("B-factor")

# 表示
plt.savefig("scatter.png")
