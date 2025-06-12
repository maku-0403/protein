import csv
import math

with open('/Volumes/pdb_res/CIF/mmCIF/00/csv_test2.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['unit','amino_number','cos_w','w'])
    
temp_list = list()
all_list = list()
new_list = list()
cos_list = list()
w_list = list()

count = 0

pass_name = "/Volumes/pdb_res/CIF/mmCIF/00/200l.cif"

with open(pass_name, "r") as f:
    for line in f:
        line = str(line)
        line = line.split(' ')
        line = [item for item in line if item != '']
        if line[0] == 'ATOM' and (line[3]=='C' or line[3]=='CA' or line[3]=='N'):
            temp_list=[line[6],line[8],line[3],line[10:13]]
            all_list.append(temp_list)
for i in range(0,len(all_list)-3):
    if all_list[i][2] == 'CA' and  all_list[i + 1][2] == 'C' and all_list[i + 2][2] == 'N' and all_list[i + 3][2] == 'CA' and abs(int(all_list[i + 1][1]) - int(all_list[i][1])) <= 1 and abs(int(all_list[i + 2][1]) - int(all_list[i + 1][1])) <= 1 and abs(int(all_list[i + 3][1]) - int(all_list[i + 2][1])) <= 1:
        
        temp = [all_list[i][0],all_list[i][1],all_list[i][3]]
        new_list.append(temp)
        
        temp = [all_list[i + 1][0],all_list[i + 1][1],all_list[i + 1][3]]
        new_list.append(temp)
        
        temp = [all_list[i + 2][0],all_list[i + 2][1],all_list[i + 2][3]]
        new_list.append(temp)
        
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
        
        cos_w = abs((aNx * bNx + aNy * bNy + aNz *bNz) / (aN * bN))
        cos_list.append(cos_w)
        
        w = math.degrees(math.acos(cos_w))
        w_list.append(w)
        
        count += 1
        
with open('/Volumes/pdb_res/CIF/mmCIF/00/csv_test2.csv', 'a') as f:
    writer = csv.writer(f)
    for i in range(0,count):
        temp_list = [new_list[i*3][0],new_list[i*3][1],cos_list[i],w_list[i]]
        writer.writerow(temp_list)
print(new_list)