import csv
import os
import glob

w_res = list()
degree = int(input('How many degree(s)? : '))
a_number = int(input('How many number(s) of amino acid? : '))

# mmCIFルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/Volumes/pdb_res/CIF/cif_to_csv/all_csv_cosw"

# 全ての.cifファイルのフルパスを再帰的に取得
csv_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

save_path_pool = ['0.5-1.0Å/', '1.0-1.5Å/', '1.5-2.0Å/', '2.0-2.5Å/', '2.5-3.0Å/', '3.0-3.5Å/', '3.5-4.0Å/', '4.0-Å/']
save_file_pool = ['0-5%', '5-10%', '10-15%', '15-20%', '20-40%', '40-60%', '60-80%', '80-100%']

path_name = '/Volumes/pdb_res/CIF/csv_to_graph_data/all_csv_resolution/all_date'

# 各PDB ID用のCSVファイルを作成
for save_path in save_path_pool:
    for save_file in save_file_pool:
        with open(f"{path_name}/each_PDBid/{save_path}{save_file}.csv", 'w') as f:
            writer = csv.writer(f)

# カウント用リストを初期化
counts = [[0 for _ in range(8)] for _ in range(8)]
process_count = 0

for csv_path in csv_files:
    # 処理中のパスを表示
    process_count += 1
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(csv_files)}")
    with open(csv_path) as f:
        reader = csv.reader(f)
        l = [row for row in reader]
        if len(l[0]) == 2:
            file_name = l[0][1]
            count = 0
            if len(l) - 6 >= a_number:
                for j in range(6, len(l)):
                    if float(l[j][3]) > degree:
                        count += 1
                w_rate = count / (len(l) - 6) * 100
                if l[1][1] != '.' and l[1][1] != 'unknown':
                    if float(l[1][1]) < 5:
                        temp = [float(l[1][1]), w_rate, file_name, l[3][1]]
                        w_res = temp
            for i, (low, high) in enumerate([(0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.0), (4.0, float('inf'))]):
                if low <= w_res[0] < high:
                    for j, (low_rate, high_rate) in enumerate([(0, 5), (5, 10), (10, 15), (15, 20), (20, 40), (40, 60), (60, 80), (80, 100)]):
                        if low_rate <= w_res[1] < high_rate:
                            counts[i][j] += 1
                            with open(f"{path_name}/each_PDBid/{save_path_pool[i]}{save_file_pool[j]}.csv", 'a') as f:
                                writer = csv.writer(f)
                                writer.writerow([w_res[2]])

# 合計を計算
totals = [sum(count) for count in counts]

# 割合を計算
def rate(a, b):
    if b == 0:
        return 0
    return f"{float(a / b * 100):.2f}%"

rates = [[rate(counts[i][j], totals[i]) for j in range(8)] for i in range(8)]

# CSVデータを保存
with open(f"{path_name}/CSV_data.csv", 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['rate', '0.5-1.0Å', '1.0-1.5Å', '1.5-2.0Å', '2.0-2.5Å', '2.5-3.0Å', '3.0-3.5Å', '3.5-4.0Å', '4.0-Å'])
    for j, save_file in enumerate(save_file_pool):
        writer.writerow([save_file] + [counts[i][j] for i in range(8)])
    for j, save_file in enumerate(save_file_pool):
        writer.writerow([save_file] + [rates[i][j] for i in range(8)])
    writer.writerow(['sum'] + totals)

print(len(w_res))