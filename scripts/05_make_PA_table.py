## python / pandas implementation for generating presence absence tables from the BLAST outputs generated previously
## depends on some of the helper functions defined in /extra/

import pandas as pd
from PythonHelperFuncs import create_df_from_blast
from PythonHelperFuncs import create_LV_dist_df

## TO DO: turn these into cmd line args
blast_output_file = "outputs/test.tsv"
annotation_file = "extra/gene_annotations.csv"
min_pident, min_qcov, verbose = (90, 85, True)
outgroup_identity_threshold = 98 ## bit of a convuluted stat, but we want to join any category that has too much similarity with the outgroup

## Importing the BLAST files
df = create_df_from_blast(blast_output_file, annotation_file, min_pident, min_qcov, verbose)
lv_df = create_LV_dist_df(blast_output_file, annotation_file, min_pident, min_qcov, target_cluster="blp", create_plots=False, verbose = verbose)

## need to drop any variants conditionally, where the 
## similarity between clusters is greater than within
## or above the threshold. The remaining variants will reflect our "real" biological variants

lv_df_kept = lv_df[(lv_df["ingroup_mean_ratio"] > lv_df["outgroup_mean_ratio"]) & (lv_df["outgroup_mean_ratio"] <= outgroup_identity_threshold)]

kept_variants = lv_df_kept["variant"].unique()

## need to create a version of df that has only variants kept in the previous step
## variants that are filtered need to be collapsed into a single variant 

## doing this by iterating over every contig. if it has multiple variants, we take the
## highest scoring that appears in the "kept_variants" list. Otherwise, we bin it into 
## "gene (non-specific)"

print(df.columns)
i = 0
for contig in df["contig"].unique():
    while i < 1:
        df_contig = df[df["contig"] == contig]
        for gene in df_contig.index.unique():
        # for row in sub_df:
            candidate_hits = df_contig[df_contig.index == gene]

            ## need to filter the variants that aren't kept
            print(candidate_hits[candidate_hits["variant"] in kept_variants])
        i += 1
