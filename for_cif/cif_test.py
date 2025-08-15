from CifFile import ReadCif

cf = ReadCif("1jsf.cif")
blk = cf.first_block()

loop = blk.GetLoop("_software.name")
cols = list(loop.keys())
rows = [list(pkt) for pkt in loop]
for i in range(0,len(rows)):
    print(rows[i]["_software.classification"])
    if rows[i]["_software.classification"] == "refinement":
        program_name = rows[i][cols.index("_software.name")]
        print('a:'+program_name)