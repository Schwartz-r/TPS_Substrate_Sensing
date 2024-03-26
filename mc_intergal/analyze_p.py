import sys
import numpy as np
import glob
import math
import pandas as pd
from read_write_pdb import PDB

def get_in():
    caver_i = int(sys.argv[1])
    rep = int(sys.argv[2])
    dim_plus=2
    p=6
    npoints=10**p
    interval=0.2
    return caver_i, rep, dim_plus, p, npoints, interval

def dist(p1,p2):
  d = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2
  d = math.sqrt(d)
  return d
  
def isin(p1,spheres):
  for s in spheres:
    d = dist(p1,s[:3])
    if d < s[3]:
      return True
  return False

def sele_p(d1,d2,interval):
    d = d1-d2
    if abs(d) < interval:
        return 0
    elif d < 0:
        return 1
    else:
        return 2

def assign(caver_i):
    if caver_i==1:
        prot="4kux"
        otrue=[ -64.076,  19.172,  10.784]
        ofalse=[ -68.424,  17.382,  13.014]
        center=[ -66.497,  22.247,   9.016]
        dim=6
    elif caver_i==2:
        prot="4okz"
        otrue=[ -69.871, 118.155,  37.167]
        ofalse=[ -65.504, 116.295,  37.424]
        center=[ -70.454, 116.181,  41.121]
        dim=6
    elif caver_i==3:
        prot="6a2c"
        otrue=[  75.776,  12.375,  60.483]
        ofalse=[  72.972,   8.460,  58.324]
        center=[  70.941,  16.594,  59.007]
        dim=6
    elif caver_i==4:
        prot="3m01"
        otrue=[ 109.189,  42.870,  50.375]
        ofalse=[ 106.257,  39.379,  49.834]
        center=[ 110.834,  40.364,  47.267]
        dim=7
    elif caver_i==5:
        prot="6wkd"
        otrue=[  78.077, -12.370, -19.580]
        ofalse=[  79.262,  -9.684, -21.557]
        center=[  81.999, -12.805, -17.260]
        dim=7
    elif caver_i==6:
        prot="5ki0"
        otrue=[  21.649,  81.570,  10.252]
        ofalse=[  24.915,  84.013,  12.773]
        center=[  23.182,  79.063,  14.387]
        dim=6
    elif caver_i==7:
        prot="3p5r"
        otrue=[  11.436,  24.388, -10.135]
        ofalse=[  10.047,  20.489,  -7.295]
        center=[  16.148,  25.061,  -7.534]
        dim=7
    elif caver_i==8:
        prot="5gue"
        otrue=[   4.159, -13.161, -13.186]
        ofalse=[   2.509, -13.006, -16.683]
        center=[   7.511, -16.354, -14.357]
        dim=6
    elif caver_i==9:
        prot="1n21"
        otrue=[  71.245,  81.700,  55.549]
        ofalse=[  66.453,  82.442,  56.787]
        center=[  69.719,  77.875,  53.674]
        dim=6
    elif caver_i==10:
        prot="2ong"
        otrue=[  24.093,  55.361, -51.842]
        ofalse=[  25.362,  50.657, -51.411]
        center=[  21.971,  56.013, -54.274]
        dim=6
    elif caver_i==11:
        prot="3v1v"
        otrue=[  21.727,  15.157,  -3.243]
        ofalse=[  22.846,  11.486,  -6.114]
        center=[  25.084,  16.210,  -0.301]
        dim=6
    elif caver_i==12:
        prot="5nx7"
        otrue=[  -2.103, -20.365,  25.941]
        ofalse=[   1.301, -18.262,  23.313]
        center=[  -4.113, -17.589,  24.594]
        dim=6
    elif caver_i==13:
        prot="5uv1"
        otrue=[  -7.931,  11.471, -21.474]
        ofalse=[ -12.584,   9.452, -21.771]
        center=[  -7.757,  13.867, -17.989]
        dim=7
    if caver_i==14:
        prot="1jfg"
        otrue=[ -20.406, 105.506,  14.004]
        ofalse=[ -22.465, 108.265,  10.886]
        center=[ -21.290, 109.185,  15.801]
        dim=7
    elif caver_i==15:
        prot="5dz2"
        otrue=[  27.422, -24.772,  24.702]
        ofalse=[  31.678, -21.766,  23.018]
        center=[  28.192, -24.455,  20.962]
        dim=7
    elif caver_i==16:
        prot="6o9p"
        otrue=[ -13.875,  10.058,  21.604]
        ofalse=[ -15.387,  15.193,  20.361]
        center=[ -14.587,   9.891,  17.137]
        dim=7
    elif caver_i==17:
        prot="6w26"
        otrue=[  -0.075,   4.917, -51.229]
        ofalse=[   4.020,   2.599, -49.297]
        center= [   2.772,   2.806, -53.695]
        dim=7
    elif caver_i==18:
        prot="7kj8"
        otrue=[   7.339,  -3.129,  20.461]
        ofalse=[  10.884,   0.766,  19.118]
        center=[   6.972,  -0.279,  20.384]
        dim=7
    elif caver_i==19:
        prot="7oc4"
        otrue=[ -19.393,   9.972, -38.187]
        ofalse=[ -15.471,   7.042, -39.300]
        center=[ -15.153,  11.597, -36.894]
        dim=7
    elif caver_i==20:
        prot="7ofl"
        otrue=[  12.399,  12.120,  49.007]
        ofalse=[  11.474,  11.234,  54.253]
        center=[  11.526,  14.250,  51.181] 
        dim=7
    elif caver_i==21:
        prot="7y9g"
        otrue=[  38.406, -40.868,   3.133]
        ofalse=[  35.940, -37.573,   5.901]
        center=[  34.015, -43.209,   3.760]
        dim=7
    elif caver_i==22:
        prot="7y88"
        otrue=[  57.068,  31.290,  52.074]
        ofalse=[  58.371,  31.915,  54.799]
        center=[  53.062,  29.520,  51.698]
        dim=7
    elif caver_i==23:
        prot="7zrn"
        otrue=[  70.216,  -3.739, -14.599]
        ofalse=[  70.369,  -2.870, -19.618]
        center=[  71.390,  -8.156, -17.157]
        dim=7
    elif caver_i==24:
        prot="8h6u"
        otrue=[  -9.976, -60.673, -70.557]
        ofalse=[ -10.828, -55.581, -70.251]
        center=[  -8.122, -58.899, -65.902]
        dim=7
    return prot,otrue,ofalse,center,dim

def get_points(caver_i,rep0,interval,dim_plus,npoints,p):
    
    fname=(f"caver_report_{caver_i}_{rep0}.csv")
    line="power,prot,i,rep,in_p1,in_p2,ratio\n"
    with open(fname,"w") as f:
        f.write(line)
    
    prot,otrue,ofalse,center,dim = assign(caver_i)
    dim+=dim_plus
    idir="../caver_output_si/"+str(caver_i)+"/data/clusters_timeless/"#tun_cl_<3-digit-int>_1.pdb
    pdbs=glob.glob(idir+"tun_cl_*_1.pdb")
    spheres=[]
    for pdbf in pdbs:
      pdb = PDB(pdbf)
      atoms = pdb.get_atoms(to_dict=False)
      for atom in atoms:
        spheres.append([atom.x,atom.y,atom.z,atom.tempFactor])
    # loose order of spheres to remove duplicates
    spheres = list(set(tuple(i) for i in spheres))
    
    low=[i-dim for i in center]
    high=[i+dim for i in center]
    size=(npoints,3)
    
    for rep in range(rep0,rep0+5):
        points = np.random.uniform(low=low, high=high, size=size)

        df = pd.DataFrame(points, columns = ['x','y','z'])
        
        df["distpa"]=df.apply(lambda p: dist([p.x, p.y, p.z],otrue), axis=1)
        df["distpb"]=df.apply(lambda p: dist([p.x, p.y, p.z],ofalse), axis=1)
        df["in"]=df.apply(lambda p: isin([p.x, p.y, p.z],spheres), axis=1)
        df["p"]=df.apply(lambda p: sele_p(p.distpa,p.distpb,interval), axis=1)
        s=df[["in","p"]].loc[(df["in"]==True) & (df["p"] > 0)].groupby("p").size()
        in_p1=s[1]
        in_p2=s[2]
        ratio=(s[1]/s[2]).round(2)
        with open(fname,"a") as f:
            line=str(p)+","+str(prot)+","+str(caver_i)+","+str(rep)+","
            line+=str(in_p1)+","+str(in_p2)+","+str(ratio)+"\n"
            f.write(line)
        df.to_csv("case"+str(caver_i)+"_rep"+str(rep)+"_"+str(npoints)+".csv")
        
def main():
    caver_i, rep0, dim_plus, p, npoints, interval = get_in()
    get_points(caver_i,rep0,interval,dim_plus,npoints,p)

if __name__ == "__main__":
    main()
