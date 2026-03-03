## This script searches all BLAST results identified in blast_outputs.txt as belonging 
## to the same gene, and then reports the number of unique sequences identified in the genome 
## ideally, you'll only have 2-3 variants per identified gene
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

blast_outputs_file = "outputs/blast_outputs.tsv"
df = pd.read_csv(blast_outputs_file, delimiter="\t")

## used for color plotting later, splitting into major and minor clusters
annotations = pd.read_csv("extra/gene_annotations.csv")
annotations[["cluster", "subcluster"]] = annotations["cluster"].str.split(" ", expand=True)

df[["GCA_ID", "qseqid"]] = df["qseqid"].astype(str).str.split("|", expand=True) ## WGS sequences
df[["gene", "sseqid"]] = df["sseqid"].str.split("|", expand=True)[[0,1]] ## reference sequences

# ## filtering data to only retain matching hits above each threshold  
min_pident = 90
min_qcov = 85

df = df.loc[~((df["pident"] >= min_pident) & (df["qcovs"] >= min_qcov))]

gene_seqs = df["sseqid"].unique()

gene = []
seq_id = []
genome = []
WGS_sequence = []
reference_sequence = []
reference_id = []

for gene_seq in gene_seqs:
    matches = df[df["sseqid"] == gene_seq]
    for index, match in matches.iterrows():
        gene.append(match.gene)
        seq_id.append(gene_seq)
        genome.append(match["GCA_ID"])
        WGS_sequence.append(match["qseq"])
        reference_id.append(match["sseqid"])
        reference_sequence.append(match["sseq"])

out_df = pd.DataFrame(data={"gene" : gene, "variant" : seq_id, "genome" : genome, "sequences" : WGS_sequence, "reference_sequences" : reference_sequence})
out_df = out_df.sample(5000)
out_df = out_df.groupby("variant", as_index=False).agg({"gene" : "first", "sequences" : "unique", "reference_sequences" : "unique"})
out_df["sequences"] = out_df["sequences"].apply(lambda x: len(x))
out_df["reference_sequences"] = out_df["reference_sequences"].apply(lambda x: len(x))
out_df["sequence_diff"] = out_df["sequences"] - out_df["reference_sequences"]

## appending annotations / clusters for coloring
out_df = out_df.set_index("gene").join(annotations.set_index("gene"), lsuffix="_left", rsuffix="")

## variant with only blp data
blp_df = out_df[out_df["cluster"] == "blp"]

sns.relplot(data=out_df, x = "sequences", y = "sequence_diff", hue="cluster", style = "functions")
sns.relplot(data=blp_df, x = "sequences", y = "sequence_diff", hue="functions")
plt.show()