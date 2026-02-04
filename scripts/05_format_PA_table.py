## Data generated during the pipeline is not stored with metadata. For later analysis, we need to 
## reattach the strain and taxon IDs to the relevant gene presence absence profile. Contigs are 
## searched using BioPython and then the relevant entry is downloaded and parsed. Realistically,
## this should only work with data uploaded and formatted like the Croucher dataset. 
## This is NOT a generalizable solution to this problem. It ONLY formats the test data generated
## by this specific pipeline on this specific dataset. Use at your own risk.

import pandas as pd
from Bio import Entrez
import os

## Entrez stuff
Entrez.email = "add_an_email@example.url"
Entrez.tool = "MightBeNiceToNameThis"

## loading the table
df_pa = pd.read_csv("outputs/BacteriocinLoci_pa_table_2026-Feb-04.csv", index_col="Unnamed: 0")

Isolate_Vec = [] ## what we actually need
Strain_Vec = [] ## not necessary explicitly, but nice for later data validation


for i in range(len(df_pa)):
    test_gca = df_pa.GCA_ID.iloc[i]

    print("Checking metadata for assembly " + str(test_gca))
    stream = Entrez.esearch(db='assembly', term = test_gca, retmax="40")
    record1 = Entrez.read(stream)

    stream = Entrez.esummary(db='assembly', id = record1['IdList'][0])
    record2 = Entrez.read(stream)
    # print(record2["DocumentSummarySet"])
    Isolate = record2["DocumentSummarySet"]["DocumentSummary"][0]["Biosource"]["Isolate"]
    Strain = record2["DocumentSummarySet"]["DocumentSummary"][0]["Biosource"]["InfraspeciesList"][1]["Sub_value"]
    
    if len(Isolate) == 0:
        Isolate = Strain ## handling unlabeled isolates the way they seem to be handled in other metadata (cf. the reference metadata file in this repo)
    
    Isolate_Vec.append(Isolate)
    Strain_Vec.append(Strain)



old_cols = list(df_pa.columns)
## tidying the columns and putting the new IDs in the front
old_cols.reverse()
old_cols.append("Strain")
old_cols.append("Isolate")
old_cols.reverse()

df_pa["Isolate"] = Isolate_Vec
df_pa["Strain"] = Strain_Vec
df_pa = df_pa[old_cols]

print("Succesfully generated labeled PA table! Saving to:\n outputs/labeled_pa_table \n This version can help you compare Strain and Isolate labels with reference data! \n")
print(df_pa)

if not os.path.exists("outputs"):
    os.makedirs("outputs")

print("Saving the table!")
df_pa.to_csv("outputs/labeled_pa_table.csv")


## Version formatted for EvoScope / PastML
print("\n\n I'm also going to make you a version formatted for PastML. You'll find it at: \n outputs/pastML_format_pa_table.csv \n This is probably what you actually want for future analysis!")
df_pastml = df_pa
df_pastml = df_pastml.rename(columns={"Isolate" : "ID"})
df_pastml = df_pastml.drop(columns=["Strain", "GCA_ID"])
df_pastml = df_pastml.set_index("ID")
df_pastml.to_csv("outputs/pastML_format_pa_table.csv")