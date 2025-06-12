import csv
import os
import glob
from pathlib import Path
import pandas as pd

# sphereルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/Volumes/pdb_res/CIF/cif_to_csv/all_csv_temperature/sphere"

# 全ての.csvファイルのフルパスを再帰的に取得
cif_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

#１つのCSVファイルごとに再帰的に繰り返し
for csv_path in cif_files:
    while True:
        path = Path(csv_path)
        pdb_id = path.stem
        
        #処理中のパスを表示
        print(f"処理中: {csv_path}")
        
        data_list = list()
        del_list = list()
        
        with open(csv_path,"r") as f:
            for line in f:
                data_list.append(line.strip().replace("'", "").split(','))

        for i in range(2,len(data_list)):
            if data_list[i][2] == "CA" and len(data_list) - i >= 6:
                j = 0
                while j < 4:
                    if float(data_list[i+j][4]) > 0.3:
                        temp_list = [data_list[i][0],data_list[i][1]]
                        del_list.append(temp_list)
                        break
                    else:
                        j += 1
                i += 3

        if len(del_list) == 0:
            break

        data_list = list()

        cos_path = "/Volumes/pdb_res/CIF/cif_to_csv/all_csv_cosw/"+pdb_id+".csv"

        with open(cos_path,"r") as f:
            for line in f:
                data_list.append(line.strip().replace("'", "").split(','))
        
        if len(data_list) == 1:
            break
        
        cos_df = pd.DataFrame(data=data_list[6:], columns=data_list[5])
        
        to_drop = pd.MultiIndex.from_tuples(del_list, names=['unit', 'amino_number'])

        cos_df = (
            cos_df
                .set_index(['unit', 'amino_number'])
                .drop(index=to_drop, errors='ignore')   # errors='ignore' で存在しない組み合わせは無視
                .reset_index()
        )


        
        out_path = "/Volumes/pdb_res/CIF/allcsv_to_retouch/0.3/"+pdb_id+".csv"
        
        with open(out_path, 'w') as f:
            writer = csv.writer(f)
            for i in range(0,6):
                writer.writerow(data_list[i])
        
        cos_df.to_csv(out_path, mode="a", header=False, index=False)
        
        break