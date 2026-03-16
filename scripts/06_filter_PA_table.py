## takes the output of step 5 and returns a PA table only keeping 
## genes / variants we want
## also filters over the FASTA files and merges duplicate sequences

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

PA_table_file = "outputs/PA_table_with_variants.csv"
annotation_file = "extra/gene_annotations.csv"
target_cluster = "blp" ## only want genes from this group
seq_similarity_threshold = .90 ## maximum percent identity for variants before we merge them


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

## checks if an identical sequence already exists, filteres it
gene_var_dict = {}
for gene_var in kept_genes:

    handle = load_gene_var_fasta(gene_var)

    for record in SeqIO.parse(handle, "fasta"):
        if str(record.id).strip().lower() == gene_var.strip().lower():
            flag = True
            for key_rec in gene_var_dict.keys():
                if record.seq == key_rec:
                    gene_var_dict[record.seq].append(gene_var)
                    flag = False
            if flag == True:
                gene_var_dict[record.seq] = [gene_var]

variants_to_keep = []
for var in gene_var_dict.keys():
    # if len(gene_var_dict[var]) > 1:
    variants_to_keep.append(sorted(gene_var_dict[var])[0])


## filters variants that are too similar, groups them together
variants_to_merge = {}
for gene_var_focal in variants_to_keep:
    
    gene_foc, var_foc = gene_var_focal.split("|")
    handle_foc = load_gene_var_fasta(gene_var_focal)
    records_foc = SeqIO.parse(handle_foc, "fasta")

    for record_foc in records_foc:
        # print(str(record_foc.id).strip().lower())
        # print(gene_var_focal.strip().lower(), "\n")
        if str(record_foc.id).strip().lower() == gene_var_focal.strip().lower():

            for gene_var_compare in variants_to_keep:
                gene_comp, var_comp = gene_var_compare.split("|")

                if gene_foc == gene_comp:

                    if var_foc != var_comp:

                        handle_comp = load_gene_var_fasta(gene_var_compare)
                        records_comp = SeqIO.parse(handle_comp, "fasta")
                        
                        for record_comp in records_comp:
                            # print(record_comp.description)
                            if str(record_comp.id).strip().lower() == gene_var_compare.strip().lower():

                                if lv.ratio(record_foc.seq, record_comp.seq) > seq_similarity_threshold:
                                    if gene_foc in variants_to_merge.keys():
                                        variants_to_merge[gene_comp].append(var_foc)
                                        variants_to_merge[gene_comp].append(var_comp)
                                    else:
                                        variants_to_merge[gene_comp] = [var_foc, var_comp]
variants_to_discard = []

print(variants_to_merge.keys())
for gene in variants_to_merge.keys():
    new_column_id = gene + "|" + "+".join(set(variants_to_merge[gene]))
    for x in set(variants_to_merge[gene]):
        variants_to_keep.remove(gene+"|"+x)
    df[new_column_id] = df[list([gene+"|"+x for x in set(variants_to_merge[gene])])].sum(axis=1).clip(0, 1)
    variants_to_keep.append(new_column_id)

df = df[sorted(list(set(variants_to_keep)))]
print(df.columns)
df = df.reset_index(names="GCA_ID")
df.to_csv("outputs/" + target_cluster + "_PA_table_with_filtered_variants_threshold_" +str(seq_similarity_threshold) +".csv")