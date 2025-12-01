import csv
import os
import glob
import pandas as pd

# CSVを読み込み（ヘッダーなし）
df = pd.read_csv("/srv/shared/pdb_resolution_list.csv", header=None)

# 1列目をキー、2列目を値として辞書に変換
mapping = pd.Series(df[1].values, index=df[0]).to_dict()

w_res = list()
while(True):
    input_setting = input("1:RUN 2:setting_params (initial:10 degrees, 100 amino number) : ")

    if input_setting == '1':
        degree = 10
        input_amino_number = 100
        break
    elif input_setting == '2':
        degree = int(input('How many degree(s)? : '))
        input_amino_number = int(input('How many number(s) of amino acid? : '))
        break

# mmCIFルートディレクトリ（あなたの環境に合わせて変更）
root_dir = input("CSV root directory: ")

# 全ての.cifファイルのフルパスを再帰的に取得
csv_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

save_path_pool = ['0.5-1.0Å', '1.0-1.5Å', '1.5-2.0Å', '2.0-2.5Å', '2.5-3.0Å', '3.0-3.5Å', '3.5-4.0Å', '4.0-4.5Å', '4.5-5.0Å', '5.0Å-']
save_file_pool = ['0-5%', '5-10%', '10-15%', '15-20%', '20-40%', '40-60%', '60-80%', '80-100%']

out_dir = input("Output directory: ")

input_program_name = input('Experiment Type： ')

todays_date = input("Today's Date: ")

os.makedirs(out_dir+"/each_PDBid", exist_ok=True)
for i in range(0,len(save_path_pool)):
    os.makedirs(out_dir+"/each_PDBid/"+save_path_pool[i], exist_ok=True)

# 各PDB ID用のCSVファイルを作成
for save_path in save_path_pool:
    for save_file in save_file_pool:
        with open(f"{out_dir}/each_PDBid/{save_path}/{save_file}.csv", 'w',newline="") as f:
            writer = csv.writer(f)

# カウント用リストを初期化
counts = [[0 for _ in range(8)] for _ in range(10)]
process_count = 0
surch_count = 0
res_error_count = 0

for csv_path in csv_files:
    process_count += 1
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(csv_files)}")
    with open(csv_path) as f:
        reader = csv.reader(f)
        l = [row for row in reader]
        if len(l) != 0:
            if (len(l[0]) == 2) and (len(l) > 8):
                #search here
                if "ELECTRON" in l[3][1]:
                    csv_program_name = ''
                    if len(l[4]) > 2:
                        for i in range(1,len(l[4])):
                            if "SHELX" in l[4][i]:
                                csv_program_name = l[4][i]
                    elif len(l[4]) == 2:
                        csv_program_name = l[4][1]
                        if l[4][1] == '?':
                            print(csv_path)
                    else:
                        continue
                    if input_program_name in csv_program_name:
                        surch_count += 1
                        file_name = l[0][1]
                        if l[1][1] == '?' or l[1][1] == '' or l[1][1] == '.' or l[1][1] == 'unknown':
                            if (file_name in mapping) and (mapping[file_name] != ''):
                                resolution = float(mapping[file_name])
                                res_error_count += 1
                            else:
                                resolution = '?'
                        else:
                            resolution = float(l[1][1])
                        multi_unit_name = [l[9][0]]
                        multi_unit_number_index = [9]
                        for i in range(10,len(l)):
                            if l[i][0] != l[i-1][0]:
                                multi_unit_name.append(l[i][0])
                                multi_unit_number_index.append(i)
                        multi_unit_number_index.append(len(l))
                        for i in range(0,len(multi_unit_name)):
                            count = 0
                            w_res = list()
                            unit_name = multi_unit_name[i]
                            amino_number = int(l[multi_unit_number_index[i+1]-1][1])
                            if  amino_number >= input_amino_number:
                                for j in range(multi_unit_number_index[i], multi_unit_number_index[i+1]):
                                    if float(l[j][5]) > degree:
                                        count += 1
                                if amino_number <= 9:
                                    break
                                w_rate = count / amino_number * 100
                                if resolution != '?':
                                    temp = [resolution, w_rate, file_name, l[3][1]]
                                    w_res = temp
                            for i, (low, high) in enumerate([(0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.0), (4.0, 4.5), (4.5, 5.0), (5.0, float('inf'))]):
                                if len(w_res) == 0:
                                    break
                                if low <= w_res[0] < high:
                                    for j, (low_rate, high_rate) in enumerate([(0, 5), (5, 10), (10, 15), (15, 20), (20, 40), (40, 60), (60, 80), (80, 100)]):
                                        if low_rate <= w_res[1] < high_rate:
                                            counts[i][j] += 1
                                            with open(f"{out_dir}/each_PDBid/{save_path_pool[i]}/{save_file_pool[j]}.csv", 'a',newline="") as f:
                                                writer = csv.writer(f)
                                                writer.writerow([w_res[2],unit_name])

# 合計を計算
totals = [sum(count) for count in counts]

# 割合を計算
def rate(a, b):
    if b == 0:
        return 0
    return f"{float(a / b * 100):.2f}%"

rates = [[rate(counts[i][j], totals[i]) for j in range(8)] for i in range(10)]

# CSVデータを保存
with open(f"{out_dir}/CSV_data_electron_"+input_program_name+'_'+todays_date+".csv", 'w',newline="") as f:
    writer = csv.writer(f)
    writer.writerow(['rate', '0.5-1.0Å', '1.0-1.5Å', '1.5-2.0Å', '2.0-2.5Å', '2.5-3.0Å', '3.0-3.5Å', '3.5-4.0Å', '4.0-4.5Å', '4.5-5.0Å', '5.0Å-'])
    for j, save_file in enumerate(save_file_pool):
        writer.writerow([save_file] + [counts[i][j] for i in range(10)])
    for j, save_file in enumerate(save_file_pool):
        writer.writerow([save_file] + [rates[i][j] for i in range(10)])
    writer.writerow(['sum'] + totals)

print(res_error_count)