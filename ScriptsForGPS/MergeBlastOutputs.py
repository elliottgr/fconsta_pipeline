import pandas as pd
import os
import glob

# df = pd.read_csv("outputs/GPS_Intermediate_BIG_pa_table_2026-Jul-13.csv")
df = pd.DataFrame()
for file in glob.glob("outputs/GPS_Intermediate_*_pa_table_2026*.csv"):
    print(file, "\n")
    df2 = pd.read_csv(file, delimiter=",")
    df2 #= df2.drop("Unnamed: 0", axis = 1 )   
    df = pd.concat([df, df2])

df = df.set_index("GCA_ID")
df.to_csv("outputs/Merged_BLAST_GPS_PA_table.csv")