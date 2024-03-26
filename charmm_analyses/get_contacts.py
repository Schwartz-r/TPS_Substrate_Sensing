from read_write_pdb import PDB
import numpy as np
import sys
import pandas as pd

sizes = {"3p5r_wwt":"c20",
         "3p5r_nwt":"c20",
         "3m01_wwt":"c15",
         "3m01_nwt":"c15",
         "6a2c_wwt":"c15",
         "6a2c_nwt":"c15",
         "5uv1_wwt":"c10",
         "5uv1_nwt":"c10",
         "2ong_wwt":"c10",
         "2ong_nwt":"c10",
         "1n21_wwt":"c10",
         "1n21_nwt":"c10",
         "6o9p_wwt":"c15",
         "6o9p_nwt":"c15",
         "2j5c_wwt":"c10",
         "2j5c_nwt":"c10",
         "5gue":"c20",
         "7y88":"c20",
         "4kux":"c15",
         "4okz":"c15",
         "6wkd":"c15",
         "3v1v":"c11",
         "5nx7":"c10",
         "1jfg":"c15",
         "5dz2":"c15",
         "6w26":"c15",
         "7kj8":"c15",
         "7oc4":"c15",
         "7ofl":"c15",
         "7y9g":"c20",
         "7zrn":"c15",
         "8h6u":"c15",
         }

def handle_input():
    prot = sys.argv[1]
    pop = sys.argv[2]

def read_f(pdb_file):
    tmpfilename = ""
    pdb = PDB(pdb_file, tmpfilename)
    atoms = pdb.get_atoms(to_dict=False)
    return atoms

def get_files(prot, pop):
    ligf = "../all_mm_best2/"+prot+"/ligand_pop-"+str(pop)+".pdb"
    envf = "../all_mm_best2/"+prot+"/envi_pop-"+str(pop)+".pdb"
    ligand = read_f(ligf)
    envi = read_f(envf)
    return ligand, envi

def get_dist(atom1, atom2):
    point1 = np.array((atom1.x, atom1.y, atom1.z))
    point2 = np.array((atom2.x, atom2.y, atom2.z))
    return np.linalg.norm(point2 - point1)

def get_sigma(atom1, atom2):
    return atom1.tempFactor + atom2.tempFactor

def exclude(prot, pop):

    po = {"1":{"p":"p1",
               "o":"o1"},
          "5":{"p":"p2",
               "o":"o5"}}
    hs = {"c10":{"h1":"h11",
                 "h2":"h12"},
          "c11":{"h1":"h11",
                 "h2":"h12"},
          "c15":{"h1":"h16",
                 "h2":"h17"},
          "c20":{"h1":"h21",
                 "h2":"h22"}}
    # "c1", "c2"
    # pairs: p-c1; o-c1; o-c2; o-h1; o-h2
    pairs = []
    pairs.append((po[str(pop)]['p'], "c1"))
    pairs.append((po[str(pop)]['o'], "c1"))
    pairs.append((po[str(pop)]['o'], "c2"))
    pairs.append((po[str(pop)]['o'], hs[sizes[prot]]['h1']))
    pairs.append((po[str(pop)]['o'], hs[sizes[prot]]['h2']))
    return pairs

def check_exclude(pairs, name1, name2):
    for pair in pairs:
        if (name2.lower() == pair[0]) and (name1.lower() == pair[1]):
            return True
    return False

def get_dists(ligand,envi):
    data = []
    columns = ["name1","resn1","resi1","segn1"]
    columns += ["name2","resn2","resi2","segn2"]
    columns += ["sigma","dist","clash089","clash075"]
    for atom1 in ligand:
        for atom2 in envi:
            dist = get_dist(atom1, atom2)
            sigma = get_sigma(atom1, atom2)
            details = [atom1.name, atom1.resName, atom1.resSeq, atom1.segID]
            details += [atom2.name, atom2.resName, atom2.resSeq, atom2.segID]
            details += [sigma, dist]
            details += [sigma*0.89 > dist]
            details += [sigma*0.75 > dist]
            data.append(details)
    df = pd.DataFrame(data, columns=columns)
    return df

def get_clashes(prot, pop):
    ligand, envi = get_files(prot, pop)
    pairs = exclude(prot, pop)
    print(pairs)
    df = get_dists(ligand,envi)
    df["exclude"] = df.apply(lambda x: check_exclude(pairs, x["name1"], x["name2"]), axis=1)
    df.to_csv(f"../results/clash_{prot}_{pop}.csv")
    nclash089 = df[~df["exclude"]]['clash089'].sum()
    nclash075 = df[~df["exclude"]]['clash075'].sum()
    nexclude = df['exclude'].sum()
    return nclash089, nclash075, nexclude

def main():
    logf = 'contacts.csv'
    with open(logf,'w') as f:
        f.write('prot,pop,Nclash089,Nclash075,Nexclude\n')
    o15 = [1,5]
    for prot in sizes:
        for pop in o15:
           nclash089, nclash075, nexclude = get_clashes(prot, pop) 
           with open(logf,'a') as f:
               line = f"{prot},{pop},{nclash089},{nclash075},{nexclude}\n"
               f.write(line)

if __name__ == "__main__":
    main()
