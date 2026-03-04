## Calculating distance between sequences in the outputs of hte BLAST step
## depends on the Levenshtein distance package because I'm too lazy to do anything more complex

import Levenshtein as lv
import numpy as np
import pandas as pd
import os

blast_outputs_file = "outputs/blast_outputs.tsv"
directory = "extra/distance_matrices"
df = pd.read_csv(blast_outputs_file, delimiter="\t")

if not os.path.exists(directory):
    os.makedirs(directory)

annotations = pd.read_csv("extra/gene_annotations.csv")
annotations[["cluster", "subcluster"]] = annotations["cluster"].str.split(" ", expand=True)

df[["GCA_ID", "qseqid"]] = df["qseqid"].astype(str).str.split("|", expand=True) ## WGS sequences
df[["gene", "sseqid"]] = df["sseqid"].str.split("|", expand=True)[[0,1]] ## reference sequences

min_pident = 90
min_qcov = 85

df = df.loc[~((df["pident"] >= min_pident) & (df["qcovs"] >= min_qcov))]

gene_seqs = df["sseqid"].unique()

gene = []
seq_id = []
genome = []
WGS_sequence = []
reference_sequence = []
reference_id = []

for gene_seq in gene_seqs:
    matches = df[df["sseqid"] == gene_seq]
    for index, match in matches.iterrows():
        gene.append(match.gene)
        seq_id.append(gene_seq)
        genome.append(match["GCA_ID"])
        WGS_sequence.append(match["qseq"])
        reference_id.append(match["sseqid"])
        reference_sequence.append(match["sseq"])

out_df = pd.DataFrame(data={"gene" : gene, "variant" : seq_id, "genome" : genome, "sequences" : WGS_sequence, "reference_sequences" : reference_sequence})
# out_df = out_df.sample(5000) 

## appending annotations / clusters for coloring
out_df = out_df.set_index("gene").join(annotations.set_index("gene"), lsuffix="_left", rsuffix="")

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