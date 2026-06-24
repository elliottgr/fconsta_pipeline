## finds and returns contigs that reported no blpA presence 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

filename = "outputs/blast_outputs.tsv"
PA_file = "outputs/labeled_pa_table_blpA_merged.csv"

df = pd.read_csv(filename, sep = "\t")

df[["contig", "contigID"]] = df["qseqid"].str.split("|", n = 1, expand = True)
df["contig"] = df["contig"].str.strip()
df[["gene", "variant"]] = df["sseqid"].str.split("|", n = 1, expand = True)
df["rel_qstart"] = df["qstart"] / df["qlen"]
df["rel_qend"] = df["qend"] / df["qlen"]
df["rel_min_qstart"] = df[["rel_qstart", "rel_qend"]].min(axis=1)
df["rel_max_qend"] = df[["rel_qstart", "rel_qend"]].max(axis=1)
print(df.head())
PA = pd.read_csv(PA_file, sep = ",")


weird = PA[(PA.blpA_non_functional == 0) & (PA.blpA_functional == 0)]
weird["GCA_ID"] = weird["GCA_ID"].str.strip()


df = df.loc[df["contig"].isin(weird.GCA_ID)]
df = df[(df["gene"] == "blpA") & (df["qcovs"] <= 75) & (df["pident"] <= 85)]


## plotting the hits by coverage and pid
# fig, ax = plt.subplots(2, 1)
sns.jointplot(data = df, x="qcovs", y="pident", kind="kde")
plt.xlabel("Query coverage")
plt.ylabel("Percent Identity")
plt.title("Comparison of blpA hits that did not meet threshold", y = .9)
# plt.show()
plt.savefig("BLAST_blpA_pid_qcov.png")
# plt.figure().close()

## plotting the hits by position on contig
plt.figure()
sns.kdeplot(data = df[["rel_min_qstart", "rel_max_qend"]].rename(columns={"rel_min_qstart":"Start position", "rel_max_qend":"End position"}), bw_adjust=.5, cut=0)
plt.xlabel("Relative position on contig")
plt.savefig("BLAST_blpA_pos_on_contig.png")
plt.show()

