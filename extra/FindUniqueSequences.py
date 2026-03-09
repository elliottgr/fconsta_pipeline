## This script searches all BLAST results identified in blast_outputs.txt as belonging 
## to the same gene, and then reports the number of unique sequences identified in the genome 
## ideally, you'll only have 2-3 variants per identified gene
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from DFfromBLAST import create_df_from_blast

min_pident, min_qcov = (90, 85)

blast_outputs_file = "outputs/blast_outputs.tsv"
annotations_file = "extra/gene_annotations.csv"
out_df = create_df_from_blast(blast_outputs_file)

annotations = pd.read_csv(annotations_file)
annotations[["cluster", "subcluster"]] = annotations["cluster"].str.split(" ", expand=True)

out_df = out_df.sample(50000)
out_df.reindex(range(len(out_df)))

out_df = out_df.groupby("variant", as_index=False).agg({"gene" : "first", "sequences" : "unique", "reference_sequences" : "unique"})
out_df["sequences"] = out_df["sequences"].apply(lambda x: len(x))
out_df["reference_sequences"] = out_df["reference_sequences"].apply(lambda x: len(x))
out_df["sequence_diff"] = out_df["sequences"] - out_df["reference_sequences"]

# ## appending annotations / clusters for coloring
out_df = out_df.set_index("gene").join(annotations.set_index("gene"), lsuffix="_left", rsuffix="")

## variant with only blp data
blp_df = out_df[out_df["cluster"] == "blp"]

sns.relplot(data=out_df, x = "sequences", y = "sequence_diff", hue="cluster", style = "functions")
sns.relplot(data=blp_df, x = "sequences", y = "sequence_diff", hue="functions")
plt.show()