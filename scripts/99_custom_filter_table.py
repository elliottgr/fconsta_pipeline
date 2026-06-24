## Variant of the filter_PA_table script that simply filters everything 
## into a single gene bin EXCEPT blpA variants, which we keep seperate 
## due to their variant functionality


import pandas as pd
from Bio import SeqIO
import Levenshtein as lv

def load_gene_var_fasta(gene_var):
    gene, variant = gene_var.split("|", 1)

    ## need to try both possible download locations
    try:
        handle = open("downloaded_data/UNIPROT_gene_data/protein_sequences/" + gene + "_protein.fasta")
    except FileNotFoundError:
        handle = open("extra/sequences/" + gene +"_" + variant + "_protein.fasta")
    return handle

PA_table_file = "outputs/BacteriocinLoci_pa_table_2026-Jun-24.csv"
annotation_file = "extra/gene_annotations.csv"
target_cluster = "blp" ## only want genes from this group
seq_similarity_threshold = .00 ## maximum percent identity for variants before we merge them


df = pd.read_csv(PA_table_file, index_col = "GCA_ID").drop(columns=["Unnamed: 0"])
an_df = pd.read_csv(annotation_file)
an_df["gene_lower"] = an_df["gene"].str.lower()
an_df = an_df.dropna(subset="cluster")
out_df_list = []
kept_genes = []

seen_gene = []

for gene_var in df.columns:
    gene, variant = gene_var.split("|", 1)
    if gene in an_df["gene_lower"].to_list():
        if an_df[an_df["gene_lower"] == gene]["cluster"].item() == target_cluster:
            kept_genes.append(gene_var)

for gene_var in kept_genes:
    gene = gene_var.split("|")[0]
    if gene != "blpa":
        gene_cols = df.filter(like=gene)
        if gene not in seen_gene:
            seen_gene.append(gene)
            df[gene] = gene_cols.sum(axis=1).clip(lower = 0, upper = 1)
    else:
        seen_gene.append(gene_var)

df = df[seen_gene]
df = df.reset_index(names="GCA_ID")
df.to_csv("outputs/" + target_cluster + "_PA_table_with_binned_variants_no_threshold.csv")