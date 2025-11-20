import csv
import os
import glob
import math
import gemmi

# mmCIFルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/srv/shared/mmCIF"
out_dir = "/srv/shared/sphere"

# 全ての.cifファイルのフルパスを再帰的に取得
cif_files = glob.glob(os.path.join(root_dir, '**', '*.cif'), recursive=True)

#１つのCIFファイルごとに再帰的に繰り返し
for cif_path in cif_files:
    process_count += 1
    #処理中のパスを表示
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(cif_files)}")

    #読み込み/計算/書き込みのためのリスト類の定義
    temp_list = list()
    sphere_list = list()

    #計算可能なペプチド結合の数の変数'count'の定義
    count = 0

    #CIFファイルの読み込み
    doc = gemmi.cif.read_file(cif_path)
    blk = doc.sole_block()

    #pdb_idの抽出
    pdb_id = blk.find_value('_entry.id')

    #C,CA,Nの座標を抜き出し、「ユニット」「アミノ酸番号」「元素」「x,y,z座標」をリスト'all_list'に追加
    table = blk.find('_atom_site.', ['label_asym_id','label_seq_id','label_atom_id','B_iso_or_equiv','group_PDB'])
    if table:  # 見つかったとき
        rows = list()
        for row in table:
            if (row[4] == "ATOM") and (((row[2] == "C") or (row[2] == "CA") or (row[2] == "N"))):
                rows.append([row[0],row[1],row[2],row[3]])
        all_list = rows

    # #CIFファイルの読み込み
    # with open(cif_path, "r") as f:
    #     #上から1行ごとに読み込み、'line'として定義
    #     for line in f:
    #         #一度すべて文字として読み込み、空白を削除
    #         line = str(line)
    #         line = line.split(' ')
    #         line = [item for item in line if item != '']
    #         #PDBidを探し、変数'pdb_id'に格納
    #         if line[0] == '_entry.id':
    #             pdb_id = line[1]
    #         #もしanisouに辿り着いていない場合（変数'anisou_TF'が"F"（偽）のとき）
    #         if (line[0] == 'ATOM')  and (len(line) > 13) and (line[3]=='C' or line[3]=='CA' or line[3]=='N') and (float(line[14]) >= 0):
    #             #anisou_listに、'ユニット','アミノ酸番号','元素名','b_factor'を追加
    #             temp_list = [line[6],line[8],line[3],line[14]]
    #             sphere_list.append(temp_list)

    #抜き出した原子の個数(all_listの要素数)分ループする
    for i in range(0,len(all_list)-3):
        #CA,C,Nの順に並んでいるもののみに条件を絞る
        if all_list[i][2] == 'CA' and  all_list[i + 1][2] == 'C' and all_list[i + 2][2] == 'N' and all_list[i + 3][2] == 'CA' and all_list[i][0] == all_list[i+1][0] and all_list[i+1][0] == all_list[i+2][0] and all_list[i+2][0] == all_list[i+3][0] and min(float(all_list[i][3]),float(all_list[i+1][3]),float(all_list[i+2][3]),float(all_list[i+3][3])) >= 0:
            sphere_list.append(all_list[i])
            sphere_list.append(all_list[i+1])
            sphere_list.append(all_list[i+2])
            sphere_list.append(all_list[i+3])

    #anisouの要素数が0以上（DNA,RNAを除くため）
    if len(sphere_list) > 0:
        #ファイル名をpdbのidにし、csvファイルを新規作成
        with open(out_dir+'/'+pdb_id+'.csv', 'w') as f:
            writer = csv.writer(f)
            #'楕円体'という情報の記述
            writer.writerow(['temperature_type','sphere'])
            #'ユニット','アミノ酸番号','元素名','B_factor','球の体積'を1行ずつ全て書き出す
            writer.writerow(['unit','amino_number','element_name','B_factor','volume'])
            for i in range(0,len(sphere_list)):
                radius = math.sqrt(float(sphere_list[i][3]) / (8 * (math.pi ** 2)))
                # 球体の体積を計算
                volume = (4/3) * math.pi * (radius ** 3)
                last_list = [sphere_list[i][0],sphere_list[i][1],sphere_list[i][2],sphere_list[i][3],volume]
                writer.writerow(last_list)
    else:
        with open(out_dir+'/'+pdb_id+'.csv', 'w') as f:
            writer = csv.writer(f)
            #'楕円体'という情報の記述
            writer.writerow(['temperature_type','sphere'])
            #'ユニット','アミノ酸番号','元素名','B_factor','球の体積'を1行ずつ全て書き出す
            writer.writerow(['unit','amino_number','element_name','B_factor','volume'])