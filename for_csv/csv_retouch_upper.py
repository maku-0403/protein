import csv
import os
import glob
from pathlib import Path
import pandas as pd
import numpy as np

# sphereルートディレクトリと出力ディレクトリ
csv_root_dir = input("CSV root directory: ")
input_out_dir = input("Out directory: ")
input_threshold = input("Pacentage of cutting: ")

out_dir = input_out_dir+"/upper_"+input_threshold
os.makedirs(out_dir, exist_ok=True)

# 全ての.csvファイルのフルパスを再帰的に取得
csv_files = glob.glob(os.path.join(csv_root_dir, '**', '*.csv'), recursive=True)

process_count = 0

# 1つのCSVファイルごとに処理（保存先に同名ファイルがなければ処理）
for csv_path in csv_files:
    file_name = os.path.basename(csv_path)
    save_path = os.path.join(out_dir, file_name)
    
    process_count += 1
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(csv_files)}")

    if not os.path.exists(save_path):
        while True:
            path = Path(csv_path)
            pdb_id = path.stem
            
            data_list = list()
            temp_list = list()
            sort_list = list()
            max_data_list = list()
            del_list = list()
            
            with open(csv_path,"r") as f:
                for line in f:
                    data_list.append(line.strip().replace("'", "").split(','))
            for i in range(2,len(data_list)):
                if data_list[i][2] == "CA" and len(data_list) - i >= 6:
                    temp_list = [data_list[i][0],data_list[i][1],float(data_list[i][4])]
                    for j in range(0,4):
                        temp_list[2] = max(temp_list[2],float(data_list[i+j][4]))
                    sort_list.append(temp_list[2])
                    max_data_list.append(temp_list)
                    i += 3
            
            if len(sort_list) == 0:
                break

            #####ここで残すパーセントを指定（75を入力したら25%切り捨てという意味）#####
            threshold = np.percentile(sort_list,100-input_threshold)
            
            for i in range(0,len(sort_list)):
                if sort_list[i] > threshold:
                    temp_list = [max_data_list[i][0],max_data_list[i][1]]
                    del_list.append(temp_list)

            data_list = list()

            cos_path = csv_root_dir+"/"+pdb_id+".csv"

            try:
                with open(cos_path,"r") as f:
                    for line in f:
                        data_list.append(line.strip().replace("'", "").split(','))
            except:
                break
            
            if len(data_list) == 1:
                break
            
            if data_list[8][0] != 'unit':
                print('error')
                break
            
            cos_df = pd.DataFrame(data=data_list[9:], columns=data_list[8])
            
            to_drop = pd.MultiIndex.from_tuples(del_list, names=['unit', 'amino_number'])

            cos_df = (
                cos_df
                    .set_index(['unit', 'amino_number'])
                    .drop(index=to_drop, errors='ignore')   # errors='ignore' で存在しない組み合わせは無視
                    .reset_index()
            )
            
            with open(save_path, 'w',newline="") as f:
                writer = csv.writer(f)
                for i in range(0,9):
                    writer.writerow(data_list[i])
            
            cos_df.to_csv(save_path, mode="a", header=False, index=False)
            
            break