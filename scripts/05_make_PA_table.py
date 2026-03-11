## python / pandas implementation for generating presence absence tables from the BLAST outputs generated previously
## depends on some of the helper functions defined in /extra/

import pandas as pd
from PythonHelperFuncs import create_df_from_blast
from PythonHelperFuncs import create_LV_dist_df

## TO DO: turn these into cmd line args
blast_output_file = "outputs/test.tsv"
annotation_file = "extra/gene_annotations.csv"
min_pident, min_qcov, verbose = (90, 85, True)
outgroup_identity_threshold = 99 ## bit of a convuluted stat, but we want to join any category that has too much similarity with the outgroup

## Importing the BLAST files
df = create_df_from_blast(blast_output_file, annotation_file, min_pident, min_qcov, verbose)
lv_df = create_LV_dist_df(blast_output_file, annotation_file, min_pident, min_qcov, target_cluster="blp", create_plots=False, verbose = verbose)

## need to drop any variants conditionally, where the 
## similarity between clusters is greater than within
## or above the threshold. The remaining variants will reflect our "real" biological variants

lv_df_kept = lv_df[(lv_df["ingroup_mean_ratio"] > lv_df["outgroup_mean_ratio"]) & (lv_df["outgroup_mean_ratio"] <= outgroup_identity_threshold)]
kept_variants = lv_df_kept["variant"].unique()

## of the remaining variants, we need to merge variants that have identical sequences
## storing as a dict where keys are the original variant ID and elements are merged IDs
# merged_variants = {}
# for variant in kept_variants:
#     sequence_row = df[df["variant"] == variant]
#     sequence = sequence_row["reference_sequences"].unique()
#     if len(sequence) > 1:
#         print(sequence_row)
#         print(sequence)

## need to create a version of df that has only variants kept in the previous step
## variants that are filtered need to be collapsed into a single variant 

## doing this by iterating over every contig. if it has multiple variants, we take the
## highest scoring that appears in the "kept_variants" list. Otherwise, we bin it 

df = df.sample(1000)
for contig in df["contig"].unique():

    df_contig = df[df["contig"] == contig]
    for gene in df_contig.index.unique():

        candidate_hits = df_contig[df_contig.index == gene]
        # print(candidate_hits)
        ## need to filter the variants that aren't kept
        # try:
        # print(candidate_hits["variant"])
        filtered_candidates = candidate_hits[candidate_hits["variant"].isin(kept_variants)]

        if len(filtered_candidates) > 1:
            max_p_id = max(filtered_candidates["p_id"])
            max_flo_cov = max(filtered_candidates["flo_cov"])
            print(filtered_candidates[(filtered_candidates["p_id"] == max_p_id) & (filtered_candidates["flo_cov"] == max_flo_cov)])
            # print("Succesfully ID'd")
            i =1
        # except KeyError:
            ## need to devise a system to handle cases where there are no variants that were kept
            ## Easiest answer is to just merge them into a conjoined label that is alphabetical
            # print("No succesful ID")
            # new_variant_label = "+".join(candidate_hits["variant"].sort_values().to_list())
            
