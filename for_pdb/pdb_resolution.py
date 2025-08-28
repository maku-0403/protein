from Bio.PDB import parse_pdb_header
import os
import glob
import csv

root_dir = "/Volumes/pdb_res/PDB/pdb"
out_dir = "/Volumes/pdb_res/PDB/"

# 全ての.entファイルのフルパスを再帰的に取得
pdb_files = glob.glob(os.path.join(root_dir, '**', '*.ent'), recursive=True)

process_count = 0
save_list = list()

#１つのPDBファイルごとに再帰的に繰り返し
for pdb_path in pdb_files:
    process_count += 1
    #処理中のパスを表示
    if process_count % 100 == 0:
        print(f"処理中: {process_count} / {len(pdb_files)}")

    with open(pdb_path) as fh:
        header = parse_pdb_header(fh)

    pdb_id = header.get("idcode")
    resolution = header.get("resolution")

    temp = [pdb_id,resolution]
    save_list.append(temp)

with open(out_dir+'pdb_resolution_list.csv','w') as f:
    writer = csv.writer(f)
    for i in range(0,len(save_list)):
        writer.writerow(save_list[i])