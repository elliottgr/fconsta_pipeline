## this file takes validation outputs (BLASTing our reference sequences against themselves) and plots the pairwise identity and coverage
## Idea from Martin G :)
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# i
## importing the file
df = pd.read_csv("outputs/validation_blast_outputs.tsv", delimiter="\t")

## plotting only genes against themselves
## need to add an extra column for grouping

df[["query_gene", "query_id"]] = df["qseqid"].str.split('_', expand=True).drop(columns=2)
df[["reference_gene", "reference_id"]] = df["sseqid"].str.split('|', expand=True).drop(columns=2)
df=df.drop(columns=["qseqid", "sseqid"])

## splitting the data into gene-gene matches and gene-gene non-matches
df["matches"] = df.query_gene == df.reference_gene

sns.jointplot(data = df, x="qcovs", y="pident", hue = "matches")
plt.xlabel("Query Coverage")
plt.ylabel("Percent Identity")
plt.title("Gene x Gene BLAST Results", y=.9)
plt.show()