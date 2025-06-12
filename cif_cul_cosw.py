import csv
import math
import os
import glob

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
    all_list = list()
    new_list = list()
    cos_list = list()
    w_list = list()
    
    final_date_TF = 'F'
    program_name_TF = 'F'

    pdb_id = ''
    resolution = ''
    final_date = ''
    exp_type = ''
    program_name = ''

    #計算可能なペプチド結合の数の変数'count'の定義
    count = 0

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
            #C,CA,Nの座標を抜き出し、「ユニット」「アミノ酸番号」「元素」「x,y,z座標」をリスト'all_list'に追加
            if len(line) > 12 and line[0] == 'ATOM' and (line[3]=='C' or line[3]=='CA' or line[3]=='N'):
                temp_list=[line[6],line[8],line[3],line[10:13]]
                all_list.append(temp_list)
            #分解能
            if line[0] == '_reflns.d_resolution_high' or line[0] == '_refine.ls_d_res_high':
                if len(line) > 2:
                    resolution = line[1]
                else:
                    resolution = 'unknown'
            #最終更新日
            if final_date_TF == 'F':
                if '_pdbx_audit_revision_history' in line[0]:
                    final_date_TF = 'T'
            elif final_date_TF == 'T':
                if line[0] == '#':
                    final_date_TF = 'F'
                elif '_pdbx_audit_revision_history' not in line[0]:
                    final_date = line[5]
            #実験手法
            if line[0] == '_exptl.method':
                if len(line) > 2:
                    for i in range(1,len(line)):
                        if line[i][-1] == "'":
                            temp = line[1:i+1]
                    if len(temp) == 1:
                        exp_type = ''
                        exp_type = temp[0].replace("'",'')
                    else:
                        exp_type = ''
                        for i in range(0,len(temp)):
                            exp_type += ' ' + temp[i]
                            exp_type = exp_type.replace("'",'')
            #プログラム
            if program_name_TF == 'F':
                if '_software.name' in line[0]:
                    if len(line) > 2:
                        program_name = line[1]
                    else:
                        program_name_TF = 'T'
            elif program_name_TF == 'T':
                if len(line) > 2 and line[1] == 'refinement':
                    program_name = line[0]
                    program_name_TF = 'F'

    if pdb_id == '':
        pdb_id = 'unknown'
    if resolution == '':
        resolution = 'unknown'
    if final_date == '':
        final_date = 'unknown'
    if exp_type == '':
        exp_type = 'unknown'
    if program_name == '':
        program_name = 'unknown'

    #抜き出した原子の個数(all_listの要素数)分ループする
    for i in range(0,len(all_list)-3):
        #CA,C,Nの順に並んでいるもののみに条件を絞る
        if all_list[i][2] == 'CA' and  all_list[i + 1][2] == 'C' and all_list[i + 2][2] == 'N' and all_list[i + 3][2] == 'CA' and abs(int(all_list[i + 1][1]) - int(all_list[i][1])) <= 1 and abs(int(all_list[i + 2][1]) - int(all_list[i + 1][1])) <= 1 and abs(int(all_list[i + 3][1]) - int(all_list[i + 2][1])) <= 1:
            #各ペプチド結合の「ユニット」「アミノ酸番号」「x,y,z座標」をリスト'new_list'に追加
            temp = [all_list[i][0],all_list[i][1]]
            new_list.append(temp)
            #前3つと後3つの原子座標から外積を算出
            ax = float(all_list[i][3][0]) - float(all_list[i + 1][3][0])
            ay = float(all_list[i][3][1]) - float(all_list[i + 1][3][1])
            az = float(all_list[i][3][2]) - float(all_list[i + 1][3][2])
            
            bx = float(all_list[i + 2][3][0]) - float(all_list[i + 3][3][0])
            by = float(all_list[i + 2][3][1]) - float(all_list[i + 3][3][1])
            bz = float(all_list[i + 2][3][2]) - float(all_list[i + 3][3][2])
            
            cx = float(all_list[i + 1][3][0]) - float(all_list[i + 2][3][0])
            cy = float(all_list[i + 1][3][1]) - float(all_list[i + 2][3][1])
            cz = float(all_list[i + 1][3][2]) - float(all_list[i + 2][3][2])
            
            aNx = ay * cz - az * cy
            aNy = ax * cz - az * cx
            aNz = ax * cy - ay * cx
            
            aN = (aNx ** 2 + aNy ** 2 + aNz ** 2) ** 0.5
            
            bNx = by * cz - bz * cy
            bNy = bx * cz - bz * cx
            bNz = bx * cy - by * cx
            
            bN = (bNx ** 2 + bNy ** 2 + bNz ** 2) ** 0.5
            #外積を元にcos_wを算出し、リスト'cos_list'に追加
            cos_w = abs((aNx * bNx + aNy * bNy + aNz *bNz) / (aN * bN))
            cos_w = max(-1.0, min(1.0, cos_w))  # 値を[-1, 1]に丸める

            cos_list.append(cos_w)
            #cos_wの値を元にw角を算出し、リスト'w_list'に追加
            w = math.degrees(math.acos(cos_w))
            w_list.append(w)
            #毎回countを1増やし、csv書き出しループ処理に使用
            count += 1
    #csvで書き出し
    with open('/Volumes/pdb_res/CIF/cif_to_csv/all_csv_cosw/'+pdb_id+'.csv', 'w') as f:
        writer = csv.writer(f)
        #DNAやRNAのように計算できないものがあった場合、"NO_DATA"と出力
        if len(new_list) == 0:
            writer.writerow(['NO_DATA'])
        #正常に計算できた場合(上記以外)、すべての「ユニット」「アミノ酸数」「cos_wの値」「w角」を出力
        else:
            writer.writerow(['PDB_id',pdb_id])
            writer.writerow(['RESOLUTION',resolution])
            writer.writerow(['FINAL_DATE',final_date])
            writer.writerow(['EXPERIMENT_TYPE',exp_type])
            writer.writerow(['PROGRAM',program_name])
            writer.writerow(['unit','amino_number','cos_w','w'])
            for i in range(0,count):
                temp_list = [new_list[i][0],new_list[i][1],cos_list[i],w_list[i]]
                writer.writerow(temp_list)