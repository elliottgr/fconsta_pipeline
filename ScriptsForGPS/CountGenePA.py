import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("outputs_from_GPS/PA_with_metadata.csv")
ann = pd.read_csv("extra/gene_annotations.csv")


## Manually adding a row for blpA func and non func
ann.loc[-1] = pd.Series({"gene" : "blpA_functional", "genbank_accession" : "na", "functions" : "Transporter", "cluster" : "blp"})
ann.loc[-2] = pd.Series({"gene" : "blpA_non_functional", "genbank_accession" : "na", "functions" : "Transporter", "cluster" : "blp"})

## creating new columns for number of genes in each annotation category

threshold = .1 ## minimum proportion of contigs returning the gene before it is present
gdf = df.groupby("GCA_ID")
annotations = {}

for annotation in set(ann.functions):
    annotations[annotation] = {}
    for gca in set(df.GCA_ID): 
        annotations[annotation][gca] = 0

## Selecting genes present in the data
for gene in ann.gene:
    if gene.lower() in df.columns or gene in df.columns:
        annotation = ann.loc[ann.gene == gene, "functions"].values[0]

        ## Iterating over each assembly
        GCAs_to_check = len(set(gdf.groups.keys()))
        i = 0

        for gca in set(gdf.groups.keys()):
            
            if i % 5000 == 0:
                print("Checking assembly # ", i+1 , " of ", GCAs_to_check, "for gene ", gene, flush=True)
            i += 1

            ## try catch loop to handle capitalized and non capitalized column names
            try:
                if gdf[gene.lower()].get_group(gca).mean() > threshold:
                    annotations[annotation][gca] += 1
            except KeyError:
                try:
                    if gdf[gene].get_group(gca).mean() > threshold:
                        annotations[annotation][gca] += 1
                except:
                    ## 
                    print("It seems like there's columns in the annotation table that are incorrectly mapping onto the PA table \n This will error appears when a gene name appeared in the column headers of both tables when the script initialized, but is not present when we actually loop over the dataframe. ")



for annotation in annotations.keys():
    df[annotation] = df.GCA_ID.map(annotations[annotation])
print(len(df))

df.to_csv("outputs_from_GPS/PA_table_with_annotation_counts.csv")