## finds and returns contigs that reported no blpA presence 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

filename = "outputs/blast_outputs.tsv"
PA_file = "outputs/labeled_pa_table_blpA_merged.csv"

df = pd.read_csv(filename, sep = "\t")
print(df.columns)
df[["contig", "contigID"]] = df["qseqid"].str.split("|", n = 1, expand = True)
df["contig"] = df["contig"].str.strip()
df[["gene", "variant"]] = df["sseqid"].str.split("|", n = 1, expand = True)
PA = pd.read_csv(PA_file, sep = ",")


weird = PA[(PA.blpA_non_functional == 0) & (PA.blpA_functional == 0)]
weird["GCA_ID"] = weird["GCA_ID"].str.strip()

# print(weird.head())

# print(df["contig"].head())

df = df.loc[df["contig"].isin(weird.GCA_ID)]
df = df[(df["gene"] == "blpA") & (df["qcovs"] <= 75) & (df["pident"] <= 85)]


## plotting the hits by coverage and pid
fig, ax = plt.subplots()
# ax.scatter(df["qcovs"], df["pident"])

sns.jointplot(data = df, x="qcovs", y="pident", kind="kde")
plt.xlabel("Query coverage")
plt.ylabel("Percent Identity")
# ax.set_label()
plt.title("Comparison of blpA hits that did not meet threshold")
plt.show()
plt.savefig("BLAST_blpA_pid_qcov.png", y = .9)


## plotting the hits by position on contig

