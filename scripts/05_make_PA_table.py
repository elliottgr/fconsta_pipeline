## python / pandas implementation for generating presence absence tables from the BLAST outputs generated previously
## depends on some of the helper functions defined in /extra/

import pandas as pd
from PythonHelperFuncs import create_df_from_blast
from PythonHelperFuncs import create_LV_dist_df
from Bio import SeqIO
import glob



## TO DO: turn these into cmd line args
blast_output_file = "outputs/test.tsv"
annotation_file = "extra/gene_annotations.csv"
main_FASTA_directory = "downloaded_data/UNIPROT_gene_data/protein_sequences/" ## needed for data validation
supplemental_FASTA_directory = "extra/sequences/"
min_pident, min_qcov, verbose = (90, 95, True)
outgroup_identity_threshold = 99 ## bit of a convuluted stat, but we want to join any category that has too much similarity with the outgroup

## Importing the BLAST files
df = create_df_from_blast(blast_output_file, annotation_file, min_pident, min_qcov, verbose)
lv_df = create_LV_dist_df(blast_output_file, annotation_file, min_pident, min_qcov, target_cluster="blp", create_plots=False, verbose = verbose, generate_new_table = False)

## need to drop any variants conditionally, where the 
## similarity between clusters is greater than within
## or above the threshold. The remaining variants will reflect our "real" biological variants

lv_df_kept = lv_df[(lv_df["ingroup_mean_ratio"] > lv_df["outgroup_mean_ratio"]) & (lv_df["outgroup_mean_ratio"] <= outgroup_identity_threshold)]
kept_variants = lv_df_kept["variant"].unique()

## need to create a version of df that has only variants kept in the previous step
## variants that are filtered need to be collapsed into a single variant 

## doing this by iterating over every contig. if it has multiple variants, we take the
## highest scoring that appears in the "kept_variants" list. Otherwise, we bin it 
i = 0
df = df.sample(1000)
for contig in df["contig"].unique():

    df_contig = df[df["contig"] == contig]
    reference_fasta_seqs = {}
    for gene in df_contig.index.unique():
        candidate_hits = df_contig[df_contig.index == gene]
        filtered_candidates = candidate_hits[candidate_hits["variant"].isin(kept_variants)]
        if len(filtered_candidates) > 0:
            max_p_id = max(filtered_candidates["p_id"])
            max_flo_cov = max(filtered_candidates["sub_cov"])
            filtered_candidates = filtered_candidates[(filtered_candidates["p_id"] == max_p_id) & (filtered_candidates["sub_cov"] == max_flo_cov)]

            ## handling case where there are multiple results with identical metrics AND reference sequences
            ## assigns a label that is a combination of the variant IDs
            if (len(filtered_candidates) > 1) & (len(filtered_candidates.reference_sequences.unique()) == 1):
                output_id = gene +"|" + "+".join(filtered_candidates["variant"].sort_values().unique())

            ## handling case where metrics are identical, but reference sequences differ
            ## in this case, assigns label of match with longest reference
            if len(filtered_candidates.reference_sequences.unique()) > 1:
                filtered_candidates["ref_len"] = filtered_candidates.reference_sequences.map(len)
                filtered_candidates = filtered_candidates[filtered_candidates["ref_len"] == max(filtered_candidates["ref_len"])]
                output_id = gene + "|" + str(filtered_candidates["variant"].astype("string").item())
                print(output_id) 

        # except KeyError:
            ## need to devise a system to handle cases where there are no variants that were kept
            ## Easiest answer is to just merge them into a conjoined label that is alphabetical
            # print("No succesful ID")
            # new_variant_label = "+".join(candidate_hits["variant"].sort_values().to_list())
            
