## This script searches all BLAST results identified in blast_outputs.txt as belonging 
## to the same gene, and then reports the number of unique sequences identified in the genome 
## ideally, you'll only have 2-3 variants per identified gene
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

blast_outputs_file = "outputs/test.tsv"
df = pd.read_csv(blast_outputs_file, delimiter="\t")


df[["GCA_ID", "qseqid"]] = df["qseqid"].str.split("|", expand=True) ## WGS sequences
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
# out_df = out_df.sample(50000)
out_df = out_df.groupby("variant", as_index=False).agg({"gene" : "first", "sequences" : "unique", "reference_sequences" : "unique"})
out_df["sequences"] = out_df["sequences"].apply(lambda x: len(x))
out_df["reference_sequences"] = out_df["reference_sequences"].apply(lambda x: len(x))

fig, ax = plt.subplots()
ax.scatter(out_df["sequences"], (out_df["sequences"] - out_df["reference_sequences"]))

ax.set_xlabel("# of observed sequences")
ax.set_ylabel("# of observed sequences - # of reference sequences")
plt.show()