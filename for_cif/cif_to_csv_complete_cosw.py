import csv
import math
import os
import glob
import gemmi

# mmCIFルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/Volumes/pdb_res/CIF/mmCIF"
out_dir = "/Volumes/pdb_res/CIF/cif_to_csv/all_csv_cosw"

# 全ての.cifファイルのフルパスを再帰的に取得
cif_files = glob.glob(os.path.join(root_dir, '**', '*.cif'), recursive=True)

process_count = 0

#１つのCIFファイルごとに再帰的に繰り返し
for cif_path in cif_files:
    process_count += 1
    #処理中のパスを表示
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(cif_files)}")
    
    #読み込み/計算/書き込みのためのリスト類の定義
    temp_list = list()
    all_list = list()
    new_list = list()
    cos_list = list()
    w_list = list()

    #計算可能なペプチド結合の数の変数'count'の定義
    count = 0

    #CIFファイルの読み込み
    doc = gemmi.cif.read_file(cif_path)
    blk = doc.sole_block()
    
    #pdb_idの抽出
    pdb_id = blk.find_value('_entry.id')
    
    #C,CA,Nの座標を抜き出し、「ユニット」「アミノ酸番号」「元素」「x,y,z座標」をリスト'all_list'に追加
    table = blk.find('_atom_site.', ['label_asym_id','label_seq_id','label_atom_id','Cartn_x', 'Cartn_y', 'Cartn_z','group_PDB'])
    if table:  # 見つかったとき
        rows = list()
        for row in table:
            if (row[6] == "ATOM") and ((row[2] == "C") or (row[2] == "CA") or (row[2] == "N")):
                rows.append([row[0],row[1],row[2],[row[3],row[4],row[5]]])
                all_list = rows
    
    #分解能の抽出
    try:
        resolution = blk.find_value("_refine.ls_d_res_high")
    except:
        resolution = "?"
    
    #最終更新日
    try:
        loop = blk.find_loop("_pdbx_audit_revision_history.revision_date")
        final_date = loop[-1]
    except:
        try:
            final_date = blk.find_value("_pdbx_audit_revision_history.revision_date")
        except:
            final_date = "?"

    #実験手法
    try:
        loop = blk.find_loop("_exptl.method")
        exp_type = loop[-1]
    except:
        try:
            exp_type = blk.find_value("_exptl.method")
        except:
            exp_type = "?"
    
    #プログラム
    try:
        table = blk.find('_software.', ['name','classification'])
        if table:  # 見つかったとき
            for row in table:
                if row[1] == "refinement":
                    program_name = row[0]
    except:
        try:
            program_name = blk.find_value("_software.name")
        except:
            program_name = "?"
    
    #R(WORK+TEST)
    try:
        loop = blk.find_loop("_refine.ls_R_factor_obs")
        r_work_test = loop[-1]
    except:
        try:
            r_work_test = blk.find_value("_refine.ls_R_factor_obs")
        except:
            r_work_test = "?"
    
    #R(WORK)
    try:
        loop = blk.find_loop("_refine.ls_R_factor_R_work")
        r_work = loop[-1]
    except:
        try:
            r_work = blk.find_value("_refine.ls_R_factor_R_work")
        except:
            r_work = "?"
    
    #R(FREE)
    try:
        loop = blk.find_loop("_refine.ls_R_factor_R_free")
        r_free = loop[-1]
    except:
        try:
            r_free = blk.find_value("_refine.ls_R_factor_R_free")
        except:
            r_free = "?"
    
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
            cos_list.append(cos_w)
            #cos_wの値を元にw角を算出し、リスト'w_list'に追加
            cos_w = max(-1.0, min(1.0, cos_w))
            w = math.degrees(math.acos(cos_w))
            w_list.append(w)
            #毎回countを1増やし、csv書き出しループ処理に使用
            count += 1
    #csvで書き出し
    with open(out_dir+'/'+pdb_id+'.csv', 'w') as f:
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
            writer.writerow(['R_VALUE(WORK+TEST)',r_work_test])
            writer.writerow(['R_VALUE(WORK)',r_work])
            writer.writerow(['R_VALUE(FREE)',r_free])
            writer.writerow(['unit','amino_number','cos_w','w'])
            for i in range(0,count):
                temp_list = [new_list[i][0],new_list[i][1],cos_list[i],w_list[i]]
                writer.writerow(temp_list)