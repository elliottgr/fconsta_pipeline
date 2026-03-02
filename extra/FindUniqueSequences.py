## This script searches all BLAST results identified in blast_outputs.txt as belonging 
## to the same gene, and then reports the number of unique sequences identified in the genome 
## ideally, you'll only have 2-3 variants per identified gene
import pandas as pd
import matplotlib.pyplot as plt

blast_outputs_file = "outputs/test.tsv"
df = pd.read_csv(blast_outputs_file, delimiter="\t")


df[["GCA_ID", "qseqid"]] = df["qseqid"].str.split("|", expand=True)
df["gene"] = df["sseqid"].str.split("|", expand=True)[0]
# ## filtering data to only retain matching hits above each threshold  
min_pident = 90
min_qcov = 85

df = df.loc[~((df["pident"] >= min_pident) & (df["qcovs"] >= min_qcov))]

gene_seqs = df["sseqid"].unique()

gene = []
seq_id = []
genome = []
sequence = []
reference_sequence = []

for gene_seq in gene_seqs:
    matches = df[df["sseqid"] == gene_seq]
    for index, match in matches.iterrows():
        gene.append(gene_seq.split("|")[0])
        seq_id.append(gene_seq.split("|")[-1])
        genome.append(match["GCA_ID"])
        sequence.append(match["sseq"])
        reference_sequence.append(match["qseq"])

out_df = pd.DataFrame(data={"gene" : gene, "seq_id" : seq_id, "genome" : genome, "sequence" : sequence, "reference_sequence" : reference_sequence})
print(out_df.head())
out_df.groupby("gene").sequence.value_counts().unstack().plot.barh()    
plt.show()