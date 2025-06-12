import csv
import matplotlib.pyplot as plt
import os
import glob

w_res = list()
degree = int(input('How many degree(s)? : '))
a_number = int(input('How many number(s) of amino acid? : '))

# mmCIFルートディレクトリ（あなたの環境に合わせて変更）
root_dir = "/Volumes/pdb_res/CIF/cif_to_csv/all_csv_cosw"

# 全ての.cifファイルのフルパスを再帰的に取得
csv_files = glob.glob(os.path.join(root_dir, '**', '*.csv'), recursive=True)

save_path_pool = ['0.5-1.0Å/','1.0-1.5Å/','1.5-2.0Å/','2.0-2.5Å/','2.5-3.0Å/','3.0-3.5Å/','3.5-4.0Å/','4.0-Å/']
save_file_pool = ['0-5%','5-10%','10-15%','15-20%','20-40%','40-60%','60-80%','80-100%']

path_name = '/Volumes/pdb_res/CIF/csv_to_graph_data/all_csv_resolution/all_date'

for i in range(0,8):
    for j in range(0,8):
        with open(path_name+'/each_PDBid/'+save_path_pool[i]+save_file_pool[j]+'.csv', 'w') as f:
            writer = csv.writer(f)

count1_0_0 = 0
count1_0_5 = 0
count1_0_10 = 0
count1_0_15 = 0
count1_20 = 0
count1_40 = 0
count1_60 = 0
count1_80 = 0

count2_0_0 = 0
count2_0_5 = 0
count2_0_10 = 0
count2_0_15 = 0
count2_20 = 0
count2_40 = 0
count2_60 = 0
count2_80 = 0

count3_0_0 = 0
count3_0_5 = 0
count3_0_10 = 0
count3_0_15 = 0
count3_20 = 0
count3_40 = 0
count3_60 = 0
count3_80 = 0

count4_0_0 = 0
count4_0_5 = 0
count4_0_10 = 0
count4_0_15 = 0
count4_20 = 0
count4_40 = 0
count4_60 = 0
count4_80 = 0

count5_0_0 = 0
count5_0_5 = 0
count5_0_10 = 0
count5_0_15 = 0
count5_20 = 0
count5_40 = 0
count5_60 = 0
count5_80 = 0

count6_0_0 = 0
count6_0_5 = 0
count6_0_10 = 0
count6_0_15 = 0
count6_20 = 0
count6_40 = 0
count6_60 = 0
count6_80 = 0

count7_0_0 = 0
count7_0_5 = 0
count7_0_10 = 0
count7_0_15 = 0
count7_20 = 0
count7_40 = 0
count7_60 = 0
count7_80 = 0

count8_0_0 = 0
count8_0_5 = 0
count8_0_10 = 0
count8_0_15 = 0
count8_20 = 0
count8_40 = 0
count8_60 = 0
count8_80 = 0

for csv_path in csv_files:
    #処理中のパスを表示
    print(f"処理中: {csv_path}")
    with open(csv_path) as f:
        reader = csv.reader(f)
        l = [row for row in reader]
        if len(l[0]) == 2:
            file_name = l[0][1]
            count = 0
            if len(l)-6 >= a_number:
                for j in range(6,len(l)):
                    if float(l[j][3]) > degree:
                        count += 1
                w_rate = count / (len(l)-6) * 100
                if l[1][1] != '.' and l[1][1] != 'unknown':
                    if float(l[1][1]) < 5:
                        temp = [float(l[1][1]),w_rate,file_name,l[3][1]]
                        w_res = temp
            if w_res[0] >= 0.5 and w_res[0] < 1.0:
                if w_res[1] >= 0 and w_res[1] < 5:
                    count1_0_0 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[0]+save_file_pool[0]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 5 and w_res[1] < 10:
                    count1_0_5 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[0]+save_file_pool[1]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 10 and w_res[1] < 15:
                    count1_0_10 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[0]+save_file_pool[2]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 15 and w_res[1] < 20:
                    count1_0_15 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[0]+save_file_pool[3]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 20 and w_res[1] < 40:
                    count1_20 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[0]+save_file_pool[4]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 40 and w_res[1] < 60:
                    count1_40 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[0]+save_file_pool[5]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 60 and w_res[1] < 80:
                    count1_60 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[0]+save_file_pool[6]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 80 and w_res[1] < 100:
                    count1_80 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[0]+save_file_pool[7]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
            if w_res[0] >= 1.0 and w_res[0] < 1.5:
                if w_res[1] >= 0 and w_res[1] < 5:
                    count2_0_0 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[1]+save_file_pool[0]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 5 and w_res[1] < 10:
                    count2_0_5 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[1]+save_file_pool[1]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 10 and w_res[1] < 15:
                    count2_0_10 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[1]+save_file_pool[2]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 15 and w_res[1] < 20:
                    count2_0_15 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[1]+save_file_pool[3]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 20 and w_res[1] < 40:
                    count2_20 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[1]+save_file_pool[4]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 40 and w_res[1] < 60:
                    count2_40 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[1]+save_file_pool[5]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 60 and w_res[1] < 80:
                    count2_60 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[1]+save_file_pool[6]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 80 and w_res[1] < 100:
                    count2_80 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[1]+save_file_pool[7]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
            if w_res[0] >= 1.5 and w_res[0] < 2.0:
                if w_res[1] >= 0 and w_res[1] < 5:
                    count3_0_0 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[2]+save_file_pool[0]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 5 and w_res[1] < 10:
                    count3_0_5 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[2]+save_file_pool[1]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 10 and w_res[1] < 15:
                    count3_0_10 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[2]+save_file_pool[2]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 15 and w_res[1] < 20:
                    count3_0_15 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[2]+save_file_pool[3]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 20 and w_res[1] < 40:
                    count3_20 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[2]+save_file_pool[4]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 40 and w_res[1] < 60:
                    count3_40 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[2]+save_file_pool[5]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 60 and w_res[1] < 80:
                    count3_60 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[2]+save_file_pool[6]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 80 and w_res[1] < 100:
                    count3_80 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[2]+save_file_pool[7]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
            if w_res[0] >= 2.0 and w_res[0] < 2.5:
                if w_res[1] >= 0 and w_res[1] < 5:
                    count4_0_0 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[3]+save_file_pool[0]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 5 and w_res[1] < 10:
                    count4_0_5 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[3]+save_file_pool[1]+'_.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 10 and w_res[1] < 15:
                    count4_0_10 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[3]+save_file_pool[2]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 15 and w_res[1] < 20:
                    count4_0_15 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[3]+save_file_pool[3]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 20 and w_res[1] < 40:
                    count4_20 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[3]+save_file_pool[4]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 40 and w_res[1] < 60:
                    count4_40 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[3]+save_file_pool[5]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 60 and w_res[1] < 80:
                    count4_60 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[3]+save_file_pool[6]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 80 and w_res[1] < 100:
                    count4_80 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[3]+save_file_pool[7]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
            if w_res[0] >= 2.5 and w_res[0] < 3.0:
                if w_res[1] >= 0 and w_res[1] < 5:
                    count5_0_0 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[4]+save_file_pool[0]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 5 and w_res[1] < 10:
                    count5_0_5 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[4]+save_file_pool[1]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 10 and w_res[1] < 15:
                    count5_0_10 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[4]+save_file_pool[2]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 15 and w_res[1] < 20:
                    count5_0_15 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[4]+save_file_pool[3]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 20 and w_res[1] < 40:
                    count5_20 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[4]+save_file_pool[4]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 40 and w_res[1] < 60:
                    count5_40 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[4]+save_file_pool[5]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 60 and w_res[1] < 80:
                    count5_60 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[4]+save_file_pool[6]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 80 and w_res[1] < 100:
                    count5_80 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[4]+save_file_pool[7]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
            if w_res[0] >= 3.0 and w_res[0] < 3.5:
                if w_res[1] >= 0 and w_res[1] < 5:
                    count6_0_0 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[5]+save_file_pool[0]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 5 and w_res[1] < 10:
                    count6_0_5 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[5]+save_file_pool[1]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 10 and w_res[1] < 15:
                    count6_0_10 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[5]+save_file_pool[2]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 15 and w_res[1] < 20:
                    count6_0_15 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[5]+save_file_pool[3]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 20 and w_res[1] < 40:
                    count6_20 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[5]+save_file_pool[4]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 40 and w_res[1] < 60:
                    count6_40 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[5]+save_file_pool[5]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 60 and w_res[1] < 80:
                    count6_60 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[5]+save_file_pool[6]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 80 and w_res[1] < 100:
                    count6_80 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[5]+save_file_pool[7]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
            if w_res[0] >= 3.5 and w_res[0] < 4.0:
                if w_res[1] >= 0 and w_res[1] < 5:
                    count7_0_0 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[6]+save_file_pool[0]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 5 and w_res[1] < 10:
                    count7_0_5 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[6]+save_file_pool[1]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 10 and w_res[1] < 15:
                    count7_0_10 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[6]+save_file_pool[2]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 15 and w_res[1] < 20:
                    count7_0_15 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[6]+save_file_pool[3]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 20 and w_res[1] < 40:
                    count7_20 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[6]+save_file_pool[4]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 40 and w_res[1] < 60:
                    count7_40 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[6]+save_file_pool[5]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 60 and w_res[1] < 80:
                    count7_60 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[6]+save_file_pool[6]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 80 and w_res[1] < 100:
                    count7_80 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[6]+save_file_pool[7]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
            if w_res[0] >= 4.0 :
                if w_res[1] >= 0 and w_res[1] < 5:
                    count8_0_0 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[7]+save_file_pool[0]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 5 and w_res[1] < 10:
                    count8_0_5 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[7]+save_file_pool[1]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 10 and w_res[1] < 15:
                    count8_0_10 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[7]+save_file_pool[2]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 15 and w_res[1] < 20:
                    count8_0_15 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[7]+save_file_pool[3]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 20 and w_res[1] < 40:
                    count8_20 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[7]+save_file_pool[4]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 40 and w_res[1] < 60:
                    count8_40 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[7]+save_file_pool[5]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 60 and w_res[1] < 80:
                    count8_60 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[7]+save_file_pool[6]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])
                if w_res[1] >= 80 and w_res[1] < 100:
                    count8_80 += 1
                    with open(path_name+'/each_PDBid/'+save_path_pool[7]+save_file_pool[7]+'.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([w_res[2]])

one = [count1_0_0,count1_0_5,count1_0_10,count1_0_15,count1_20,count1_40,count1_60,count1_80,count1_0_0+count1_0_5+count1_0_10+count1_0_15+count1_20+count1_40+count1_60+count1_80]
two = [count2_0_0,count2_0_5,count2_0_10,count2_0_15,count2_20,count2_40,count2_60,count2_80,count2_0_0+count2_0_5+count2_0_10+count2_0_15+count2_20+count2_40+count2_60+count2_80]
three = [count3_0_0,count3_0_5,count3_0_10,count3_0_15,count3_20,count3_40,count3_60,count3_80,count3_0_0+count3_0_5+count3_0_10+count3_0_15+count3_20+count3_40+count3_60+count3_80]
four = [count4_0_0,count4_0_5,count4_0_10,count4_0_15,count4_20,count4_40,count4_60,count4_80,count4_0_0+count4_0_5+count4_0_10+count4_0_15+count4_20+count4_40+count4_60+count4_80]
five = [count5_0_0,count5_0_5,count5_0_10,count5_0_15,count5_20,count5_40,count5_60,count5_80,count5_0_0+count5_0_5+count5_0_10+count5_0_15+count5_20+count5_40+count5_60+count5_80]
six = [count6_0_0,count6_0_5,count6_0_10,count6_0_15,count6_20,count6_40,count6_60,count6_80,count6_0_0+count6_0_5+count6_0_10+count6_0_15+count6_20+count6_40+count6_60+count6_80]
seven = [count7_0_0,count7_0_5,count7_0_10,count7_0_15,count7_20,count7_40,count7_60,count7_80,count7_0_0+count7_0_5+count7_0_10+count7_0_15+count7_20+count7_40+count7_60+count7_80]
eight = [count8_0_0,count8_0_5,count8_0_10,count8_0_15,count8_20,count8_40,count8_60,count8_80,count8_0_0+count8_0_5+count8_0_10+count8_0_15+count8_20+count8_40+count8_60+count8_80]

one_count = sum(one[0:8])
two_count = sum(two[0:8])
three_count = sum(three[0:8])
four_count = sum(four[0:8])
five_count = sum(five[0:8])
six_count = sum(six[0:8])
seven_count = sum(seven[0:8])
eight_count = sum(eight[0:8])

def rate(a,b):
    if b == 0:
        return 0
    rates = float(a/b*100)
    return str(rates)+'%'

one_all = [rate(one[0],one[8]),rate(one[1],one[8]),rate(one[2],one[8]),rate(one[3],one[8]),rate(one[4],one[8]),rate(one[5],one[8]),rate(one[6],one[8]),rate(one[7],one[8])]
two_all = [rate(two[0],two[8]),rate(two[1],two[8]),rate(two[2],two[8]),rate(two[3],two[8]),rate(two[4],two[8]),rate(two[5],two[8]),rate(two[6],two[8]),rate(two[7],two[8])]
three_all = [rate(three[0],three[8]),rate(three[1],three[8]),rate(three[2],three[8]),rate(three[3],three[8]),rate(three[4],three[8]),rate(three[5],three[8]),rate(three[6],three[8]),rate(three[7],three[8])]
four_all = [rate(four[0],four[8]),rate(four[1],four[8]),rate(four[2],four[8]),rate(four[3],four[8]),rate(three[4],three[8]),rate(three[5],three[8]),rate(three[6],three[8]),rate(three[7],three[8])]
five_all = [rate(five[0],five[8]),rate(five[1],five[8]),rate(five[2],five[8]),rate(five[3],five[8]),rate(five[4],five[8]),rate(five[5],five[8]),rate(five[6],five[8]),rate(five[7],five[8])]
six_all = [rate(six[0],six[8]),rate(six[1],six[8]),rate(six[2],six[8]),rate(six[3],six[8]),rate(five[4],five[8]),rate(five[5],five[8]),rate(five[6],five[8]),rate(five[7],five[8])]
seven_all = [rate(seven[0],seven[8]),rate(seven[1],seven[8]),rate(seven[2],seven[8]),rate(seven[3],seven[8]),rate(five[4],five[8]),rate(five[5],five[8]),rate(five[6],five[8]),rate(five[7],five[8])]
eight_all = [rate(eight[0],eight[8]),rate(eight[1],eight[8]),rate(eight[2],eight[8]),rate(eight[3],eight[8]),rate(five[4],five[8]),rate(five[5],five[8]),rate(five[6],five[8]),rate(five[7],five[8])]

with open(path_name+'/CSV_data.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['rate','0.5-1.0Å','1.0-1.5Å','1.5-2.0Å','2.0-2.5Å','2.5-3.0Å','3.0-3.5Å','3.5-4.0Å','4.0-Å'])
    writer.writerow(['0-5%',count1_0_0,count2_0_0,count3_0_0,count4_0_0,count5_0_0,count6_0_0,count7_0_0,count8_0_0])
    writer.writerow(['5-10%',count1_0_5,count2_0_5,count3_0_5,count4_0_5,count5_0_5,count6_0_5,count7_0_5,count8_0_5])
    writer.writerow(['10-15%',count1_0_10,count2_0_10,count3_0_10,count4_0_10,count5_0_10,count6_0_10,count7_0_10,count8_0_10])
    writer.writerow(['15-20%',count1_0_15,count2_0_15,count3_0_15,count4_0_15,count5_0_15,count6_0_15,count7_0_15,count8_0_15])
    writer.writerow(['20-40%',count1_20,count2_20,count3_20,count4_20,count5_20,count6_20,count7_20,count8_20])
    writer.writerow(['40-60%',count1_40,count2_40,count3_40,count4_40,count5_40,count6_40,count7_40,count8_40])
    writer.writerow(['60-80%',count1_60,count2_60,count3_60,count4_60,count5_60,count6_60,count7_60,count8_60])
    writer.writerow(['80-100%',count1_80,count2_80,count3_80,count4_80,count5_80,count6_80,count7_80,count8_80])
    writer.writerow(['0-5%',rate(one[0],one[8]),rate(two[0],two[8]),rate(three[0],three[8]),rate(four[0],four[8]),rate(five[0],five[8]),rate(six[0],six[8]),rate(seven[0],seven[8]),rate(eight[0],eight[8])])
    writer.writerow(['5-10%',rate(one[1],one[8]),rate(two[1],two[8]),rate(three[1],three[8]),rate(four[1],four[8]),rate(five[1],five[8]),rate(six[1],six[8]),rate(seven[1],seven[8]),rate(eight[1],eight[8])])
    writer.writerow(['10-15%',rate(one[2],one[8]),rate(two[2],two[8]),rate(three[2],three[8]),rate(four[2],four[8]),rate(five[2],five[8]),rate(six[2],six[8]),rate(seven[2],seven[8]),rate(eight[2],eight[8])])
    writer.writerow(['15-20%',rate(one[3],one[8]),rate(two[3],two[8]),rate(three[3],three[8]),rate(four[3],four[8]),rate(five[3],five[8]),rate(six[3],six[8]),rate(seven[3],seven[8]),rate(eight[3],eight[8])])
    writer.writerow(['20-40%',rate(one[4],one[8]),rate(two[4],two[8]),rate(three[4],three[8]),rate(four[4],four[8]),rate(five[4],five[8]),rate(six[4],six[8]),rate(seven[4],seven[8]),rate(eight[4],eight[8])])
    writer.writerow(['40-60%',rate(one[5],one[8]),rate(two[5],two[8]),rate(three[5],three[8]),rate(four[5],four[8]),rate(five[5],five[8]),rate(six[5],six[8]),rate(seven[5],seven[8]),rate(eight[5],eight[8])])
    writer.writerow(['60-80%',rate(one[6],one[8]),rate(two[6],two[8]),rate(three[6],three[8]),rate(four[6],four[8]),rate(five[6],five[8]),rate(six[6],six[8]),rate(seven[6],seven[8]),rate(eight[6],eight[8])])
    writer.writerow(['80-100%',rate(one[7],one[8]),rate(two[7],two[8]),rate(three[7],three[8]),rate(four[7],four[8]),rate(five[7],five[8]),rate(six[7],six[8]),rate(seven[7],seven[8]),rate(eight[7],eight[8])])
    writer.writerow(['sum',one_count,two_count,three_count,four_count,five_count,six_count,seven_count,eight_count])

print(len(w_res))