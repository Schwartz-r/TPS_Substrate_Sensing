#!/home/qnt/renanak/miniconda3/envs/my-rdkit-env/bin/python

with open("clashes_075.txt") as f:
    lines = f.readlines()

newlines = ""
prot = ""
for l in lines:
    newline = ""
    if ("o1a" in l) or ("o2a" in l):
        prot = l.replace("\n","")
    elif "Name" in l:
        prot = ""
    else:
        newl = l.replace("(","").replace(")",",").replace(" ","")
        newline = prot + "," + newl
    newlines += newline

with open("clashes_075.csv","w") as f:
    f.write(newlines)
        
