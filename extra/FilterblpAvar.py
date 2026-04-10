## this file opens the blast outputs and determines whether a hit contained the AA substitution that renders the blpA gene non-functional
## since the blpA functional / non-functional mutation is so small, BLAST cannot reasonably filter the variants
## We'll just use the documented substitution to categorize contigs into one of two categories
## finally, the script saves the GCA_IDs so they can be joined to the presence absence table
import pandas as pd

filename = "outputs/blast_outputs.tsv"
old_PA_table = "outputs/DataForCurnagl/labeled_pa_table_binned_genes.csv"
df = pd.read_csv(filename, sep="\t")

non_functionals = []
functionals = []

for idx, row in df.iterrows():
    if "blpA" in row.stitle:
        if "GLLSKLP" in row.qseq:
            non_functionals.append(row.qseqid.split("|")[0])
        elif "GLLSFLPL" in row.qseq:
            functionals.append(row.qseqid.split("|")[0])

PA = pd.read_csv(old_PA_table, sep="\t")
blpA_func = []
blpA_non_func = []

for idx, row in PA.iterrows():
    if row.GCA_ID in functionals:
        blpA_func.append(1)
    else:
        blpA_func.append(0)

    if row.GCA_ID in non_functionals:
        blpA_non_func.append(1)
    else:
        blpA_non_func.append(0)

new_PA = PA.copy()
new_PA.insert(4, "blpA_non_functional", blpA_non_func)
new_PA.insert(5, "blpA_functional", blpA_func)
df_pastml = new_PA.copy()
df_pastml = df_pastml.rename(columns={"Isolate" : "ID"})
df_pastml = df_pastml.drop(columns=["Strain", "GCA_ID", "Unnamed: 0"])
df_pastml = df_pastml.set_index("ID")
df_pastml.to_csv("outputs/DataForCurnagl/labeled_pa_table_binned_genes_blpA_cleaned.csv")
print(new_PA.head())