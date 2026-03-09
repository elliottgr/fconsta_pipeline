## Calculating distance between sequences in the outputs of hte BLAST step
## depends on the Levenshtein distance package because I'm too lazy to do anything more complex
## we want to generate two tables per gene: one that compares all sequences of a given variant to themselves
## and another that compares all variants 

import Levenshtein as lv
import numpy as np
import pandas as pd
import os
from DFfromBLAST import create_df_from_blast

directory = "extra/distance_matrices"
if not os.path.exists(directory):
    os.makedirs(directory)

min_pident = 90
min_qcov = 85

blast_outputs_file = "outputs/blast_out.tsv"
annotations_file = "extra/gene_annotations.csv"
out_df = create_df_from_blast(blast_outputs_file, annotations_file, min_pident, min_qcov)

## variant with only blp data
blp_df = out_df[out_df["cluster"] == "blp"].groupby("gene")

## iterating over the grouped dataframe, only focusing on blp
## each sequence will have it's LS distance calc'd then saved to 
## doing a lazy df construction

gene_variant_id = {}
foc_seqs = []
tar_seqs = []
i, j = (0,0)
gvid = 0
sequences = len(blp_df["sequences"].unique())
output_matrix = np.empty((sequences, sequences), float)


for name, group in blp_df:
    print("Now calculating for gene: "+ str(name))
    file = name + "_distance_matrix.txt"
    file_no_genomes = name + "_distance_matrix_no_gene_names.txt"
    genome_dict = {}
    i = 0
    sequences = len(group["sequences"].unique())
    output_matrix = np.empty((sequences, sequences), float)

    for seq_foc in group["sequences"].unique(): ## focal sequence
        genome_dict[i] = group[group["sequences"] == seq_foc]["genome"]
        j = 0
        for seq_tar in group["sequences"].unique(): ## target sequence
            output_matrix[i, j] = lv.distance(seq_foc, seq_tar)
            j += 1
        i += 1

    with open(directory + "/" + file, "w") as f:
        f.write('\t' + str(sequences) + "\n")
        for line in range(sequences):
            line_to_write = "+".join(genome_dict[line].unique()) + "\t" + "\t".join([str(x) for x in output_matrix[line]]) + "\n"
            f.write(line_to_write)
        f.close()
    with open(directory + "/" + file_no_genomes, "w") as f2:
        f2.write('\t' + str(sequences) + "\n")
        for line in range(sequences):
            line_to_write = "sequence_" + str(line) + "\t" + "\t".join([str(x) for x in output_matrix[line]]) + "\n"
            f2.write(line_to_write)
        f2.close()
        

print("Done!!!")