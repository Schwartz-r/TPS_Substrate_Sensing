import pandas as pd
import glob

csvs = glob.glob('caver_report_*_*.csv')
dfs = []
for f in csvs:
    df = pd.read_csv(f)
    dfs.append(df)
df = pd.concat(dfs, ignore_index=True)
grouped_stats=df.groupby("i")["ratio"].agg(["mean", "std"])
grouped_stats.to_csv("caver_final.csv")
