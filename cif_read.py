import csv
with open('/Volumes/pdb_res/CIF/mmCIF/00/csv_test1.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['unit','amino_number','element_name','x','y','z'])
temp_list = list()
path_name = "/Volumes/pdb_res/CIF/mmCIF/00/200l.cif"
with open(path_name, "r") as f:
    for line in f:
        line = str(line)
        line = line.split(' ')
        line = [item for item in line if item != '']
        print(line)
        if line[0] == 'ATOM' and (line[3]=='C' or line[3]=='CA' or line[3]=='N'):
            temp_list=[line[6],line[8],line[3],line[10],line[11],line[12]]
            with open('/Volumes/pdb_res/CIF/mmCIF/00/csv_test1.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow(temp_list)
print(temp_list)