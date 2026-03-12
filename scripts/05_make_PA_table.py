## python / pandas implementation for generating presence absence tables from the BLAST outputs generated previously
## depends on some of the helper functions defined in /extra/

import pandas as pd
from PythonHelperFuncs import create_df_from_blast
from PythonHelperFuncs import create_LV_dist_df


## TO DO: turn these into cmd line args
blast_output_file = "outputs/test.tsv"
annotation_file = "extra/gene_annotations.csv"
main_FASTA_directory = "downloaded_data/UNIPROT_gene_data/protein_sequences/" ## needed for data validation
supplemental_FASTA_directory = "extra/sequences/"
min_pident, min_qcov, verbose = (90, 85, True)
outgroup_identity_threshold = 98 ## bit of a convuluted stat, but we want to join any category that has too much similarity with the outgroup

## Importing the BLAST files
df = create_df_from_blast(blast_output_file, annotation_file, min_pident, min_qcov, verbose)
lv_df = create_LV_dist_df(blast_output_file, annotation_file, min_pident, min_qcov, target_cluster="blp", create_plots=False, verbose = verbose, generate_new_table = False)
print(df.index.unique())
## some basic pre-processing, removing \t and \n characters from sequences
df["reference_sequences"] = df["reference_sequences"].str.strip().replace("\t", "").replace("\n", "").replace(" ", "")

## need to drop any variants conditionally, where the 
## similarity between clusters is greater than within
## or above the threshold. The remaining variants will reflect our "real" biological variants

lv_df_kept = lv_df[(lv_df["ingroup_mean_ratio"] >= lv_df["outgroup_mean_ratio"])]# & (lv_df["outgroup_mean_ratio"] <= outgroup_identity_threshold)]
lv_df_kept = lv_df
kept_variants = lv_df_kept["variant"].unique().tolist()
print(kept_variants)

## need to create a version of df that has only variants kept in the previous step
## variants with identical sequences need to be collapsed
print("Checking for identical reference sequences!")
variant_merge_dict = {}
for variant in kept_variants:
    ref_sequence = df[df["variant"] == variant]["reference_sequences"].unique()
    if len(df[df["reference_sequences"].isin(ref_sequence)]["variant"].unique()) > 1:
        dupe_variants = df[df["reference_sequences"].isin(ref_sequence)]["variant"].sort_values().unique()
        for dupe_variant in dupe_variants:
            merged_label = "+".join(dupe_variants)
            variant_merge_dict[dupe_variant] = merged_label
            kept_variants.append(merged_label)

# df["variant"] = df["variant"].replace(variant_merge_dict)

## doing this by iterating over every contig. if it has multiple variants, we take the
## highest scoring that appears in the "kept_variants" list. Otherwise, we bin it 
## saving to a dict that becomes a dataframe. Each key is an isolate (NOT CONTIG)

output_dict = {}

print("Generating presence / absence table!")
for contig in df["contig"].unique():

    df_contig = df[df["contig"] == contig]
    strainID = df_contig["genome"].values[0]

    for gene in df_contig.index.unique():
        candidate_hits = df_contig[df_contig.index == gene]
        filtered_candidates = candidate_hits #[candidate_hits["variant"].isin(kept_variants)]
        if len(filtered_candidates) > 0:
            max_p_id = max(filtered_candidates["p_id"])
            max_flo_cov = max(filtered_candidates["sub_cov"])
            filtered_candidates = filtered_candidates[(filtered_candidates["p_id"] == max_p_id) & (filtered_candidates["sub_cov"] == max_flo_cov)]

            ## handling case where there are multiple results with identical metrics AND reference sequences
            ## assigns a label that is a combination of the variant IDs
            for name, candidate in filtered_candidates.iterrows():
                output_gene_label = gene + "|" + candidate["variant"]
                    # else: ## Some sequences seem to break this somehow. As an arbitrary work around, I'm simply declaring it ambiguous
                    #     output_gene_label = gene + "|ambiguous"
                # else:
                    # output_gene_label = gene + "|" + filtered_candidates["variant"].item()
                
                if strainID not in output_dict.keys():
                    output_dict[strainID] = {output_gene_label : 1}
                elif output_gene_label not in output_dict[strainID]:
                    output_dict[strainID] = {output_gene_label : 1}
                else:
                    output_dict[strainID][output_gene_label] += 1
                    
output_df = pd.DataFrame.from_dict(output_dict).transpose().fillna(0).astype(int)
print("Done! Saving to outputs/Python_PA_table.csv :)")
output_df.to_csv("outputs/Python_PA_table.csv")