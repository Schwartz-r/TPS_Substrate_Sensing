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

def main():
    o15 = [1,5]
    for prot in sizes:
        for pop in o15:
            if pop == 1:
                popo = 'o2a'
            elif pop == 5:
                popo = 'o1a'
            df = pd.read_csv(f"../results/clash_{prot}_{pop}.csv")
            ress = df[(~df["exclude"]) & (df['clash089'])][["resn2","resi2","segn2"]].apply(lambda x: tuple(x), axis=1).value_counts()
            print(prot, popo)
            print(ress, '\n')

if __name__ == "__main__":
    main()
