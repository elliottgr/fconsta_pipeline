## takes the output of step 5 and returns a PA table only keeping 
## genes / variants we want
## also filters over the FASTA files and merges duplicate sequences

import pandas as pd
from Bio import Seq.IO

PA_table_file = "outputs/PA_table_with_variants.csv"
annotation_file = "extra/gene_annotations.csv"
target_cluster = "blp" ## only want genes from this group

df = pd.read_csv(PA_table_file, index_col = "GCA_ID")
an_df = pd.read_csv(annotation_file)
an_df["gene_lower"] = an_df["gene"].str.lower()
an_df = an_df.dropna(subset="cluster")

out_df_list = []
kept_genes = []

for gene_var in df.columns:
    gene, variant = gene_var.split("|", 1)
    if gene in an_df["gene_lower"].to_list():
        if an_df[an_df["gene_lower"] == gene]["cluster"].item() == target_cluster:
            kept_genes.append(gene_var)

df.to_csv("outputs/" + target_cluster + "_PA_table_with_variants.csv")