import csv
import os
import glob
import math

# mmCIFルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/Volumes/pdb_res/CIF/mmCIF"

# 全ての.cifファイルのフルパスを再帰的に取得
cif_files = glob.glob(os.path.join(root_dir, '**', '*.cif'), recursive=True)

#１つのCIFファイルごとに再帰的に繰り返し
for cif_path in cif_files:
    #処理中のパスを表示
    print(f"処理中: {cif_path}")
    
    #読み込み/計算/書き込みのためのリスト類の定義
    temp_list = list()
    anisou_list = list()

    #計算可能なペプチド結合の数の変数'count'の定義
    count = 0

    #anisouの部分の処理状態の変数'anisou_TF'と、データ内にanisouがある/計算済を証明する'anisou_comp'（初期状態はどちらも偽）
    anisou_TF = 'F'
    anisou_comp = False

    #CIFファイルの読み込み
    with open(cif_path, "r") as f:
        #上から1行ごとに読み込み、'line'として定義
        for line in f:
            #一度すべて文字として読み込み、空白を削除
            line = str(line)
            line = line.split(' ')
            line = [item for item in line if item != '']
            #PDBidを探し、変数'pdb_id'に格納
            if line[0] == '_entry.id':
                pdb_id = line[1]
            #もしanisouに辿り着いていない場合（変数'anisou_TF'が"F"（偽）のとき）
            if anisou_TF == 'F':
                #もし'_atom_site_anisotrop'の記述を見つけたら、変数'anisou_TF'を"T"（真）にする
                if '_atom_site_anisotrop' in line[0]:
                    anisou_TF = 'T'
            #もしanisouに辿り着いている場合（変数'anisou_TF'が"T"（真）のとき）
            elif anisou_TF == 'T':
                #もし'＃'（別の新しい記述ループが始まるため）の記述を見つけたら、変数'anisou_TF'を"F"（偽）にする
                if line[0] == '#':
                    anisou_TF = 'F'
                    anisou_comp = True
                #もし'_atom_site_anisotrop'が含まれていなく、元素名が'C','CA','N'で、アミノ酸番号の部分が空('.')でない場合
                elif ('_atom_site_anisotrop' not in line[0]) and (line[2]=='C' or line[2]=='CA' or line[2]=='N') and (line[6] != '.') and (float(line[8]) >= 0) and (float(line[9]) >= 0) and (float(line[10]) >= 0):
                    #anisou_listに、'ユニット','アミノ酸番号','元素名','Uの分散値'を追加
                    temp_list = [line[5],line[6],line[2],line[8:11]]
                    anisou_list.append(temp_list)

    #anisouの要素数が0以上（DNA,RNAを除くため）で、anisouの記述がある場合
    if len(anisou_list) > 0 and anisou_comp:
        #ファイル名をpdbのidにし、csvファイルを新規作成
        with open('/Volumes/pdb_res/CIF/cif_to_csv/all_csv_temperature/ellipsoid/'+pdb_id+'.csv', 'w') as f:
            writer = csv.writer(f)
            #'楕円体'という情報の記述
            writer.writerow(['temperature_type','ellipsoid'])
            #'ユニット','アミノ酸番号','元素名','U11(x軸)','U22(y軸)','U33(z軸)'を1行ずつ全て書き出す
            writer.writerow(['unit','amino_number','element_name','x-radius','y-radius','z-radius',"volume"])
            for i in range(0,len(anisou_list)):
                temp_list = [anisou_list[i][0],anisou_list[i][1],anisou_list[i][2],anisou_list[i][3][0],anisou_list[i][3][1],anisou_list[i][3][2]]
                # 各軸方向の半径を計算
                rx = math.sqrt(float(temp_list[3]))
                ry = math.sqrt(float(temp_list[4]))
                rz = math.sqrt(float(temp_list[5]))
                # 楕円体の体積を計算
                volume = (4/3) * math.pi * rx * ry * rz
                last_list = [anisou_list[i][0],anisou_list[i][1],anisou_list[i][2],anisou_list[i][3][0],anisou_list[i][3][1],anisou_list[i][3][2],volume]
                writer.writerow(last_list)