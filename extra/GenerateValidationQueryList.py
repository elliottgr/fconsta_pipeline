## For data validation, we need to generate a list of individual sequences to blast against the DB
## run this from the parent directory 
import os

sequence_dir = "extra/sequences"
query_list = open("extra/queries.txt", "w")

## Uniprot genes
with open("downloaded_data/UNIPROT_gene_data/protein_sequences/protein_sequences.fasta", "r") as uniprot_seqs:

    ## setting a dummy file for the loop
    current_file = open("downloaded_data/UNIPROT_gene_data/protein_sequences/protein_sequences.fasta", "r")

    for line in uniprot_seqs:
        if line[0] == ">":
            current_file.close()
            current_protein = line.rstrip().split(">")[1].replace("|", "_")
            print("Found protein: " + current_protein)
            current_filename = current_protein + "_protein.fasta"
            current_file = open("extra/sequences/" + current_filename, "w")
            current_file.write(">" + current_protein + "\n")
            query_list.write("extra/sequences/" + current_filename + "\n")
        else:
            current_file.write(line + "\n")
        