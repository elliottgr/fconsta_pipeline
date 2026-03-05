## I'm curious about the maximum sequence divergence between blpA isolates

import numpy as np
import pandas as pd

blast_file = "outputs/blast_outputs.tsv"
df = pd.read_csv(blast_file, delimiter="\t")

df[["GCA_ID", "qseqid"]] = df["qseqid"].astype(str).str.split("|", expand=True) ## WGS sequences
df[["gene", "sseqid"]] = df["sseqid"].str.split("|", expand=True)[[0,1]] ## reference sequences

min_pident = 90
min_qcov = 85

df = df.loc[~((df["pident"] >= min_pident) & (df["qcovs"] >= min_qcov))]
blpA_df = df[df["gene"] == "blpA"]
blpB_df = df[df["gene"] == "blpB"]
blpM_df = df[df["gene"] == "blpM"]
blpN_df = df[df["gene"] == "blpN"]

## should be several thousand
print(blpA_df)

## should be much less than several thousand
print(blpB_df)
print(blpM_df)

## looking at the number of hits per gene per genome
## trying to see if blpA appears multiple times per genome
i = 0
for GCA in df["GCA_ID"].unique():
    if i <= 5:
        ## only looking at combinations of gene, sseqid, and strain
        print(df[["gene", "sseqid"]].value_counts()["blpA"])
        print(df[["gene", "sseqid"]].value_counts()["blpB"])
        print(df[["gene", "sseqid"]].value_counts()["blpM"])
        print(df[["gene", "sseqid"]].value_counts()["blpN"])
        i+= 1