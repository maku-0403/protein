import pandas as pd

# CSVを読み込み（ヘッダーなし）
df = pd.read_csv("/Volumes/pdb_res/PDB/pdb_resolution_list.csv", header=None)

# 1列目をキー、2列目を値として辞書に変換
mapping = pd.Series(df[1].values, index=df[0]).to_dict()

print(mapping["100D"])
