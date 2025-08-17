# from CifFile import ReadCif

# cf = ReadCif("1jsf.cif")
# blk = cf.first_block()

# loop = blk.GetLoop("_software.name")
# cols = list(loop.keys())
# rows = [list(pkt) for pkt in loop]
# for i in range(0,len(rows)):
#     print(rows[i]["_software.classification"])
#     if rows[i]["_software.classification"] == "refinement":
#         program_name = rows[i][cols.index("_software.name")]
#         print('a:'+program_name)


import gemmi

# CIFファイル読み込み
doc = gemmi.cif.read_file("/Users/kuniimahan/Downloads/1jsf.cif")

# 最初のブロックを取得
blk = doc.sole_block()

# 単独
entry_id = blk.find_value('_entry.id')
print(entry_id)

# Loop
loop = blk.find_loop('_software.classification')
for row in loop:
    print(row)  # 各rowはタプル状で列データが入っている


rows = list()

table = blk.find('_atom_site.', ['label_asym_id','label_seq_id','label_atom_id','Cartn_x', 'Cartn_y', 'Cartn_z'])
if table:  # 見つかったとき
    for row in table:
        rows.append([row[0],row[1],row[2],[row[3],row[4],row[5]]])
print(rows)

# ./KIOXIA/mmCIF