## Helper script to parse the BLAST outputs and turn them into a nice dataframe

import pandas as pd

def create_df_from_blast(blast_outputs_file, annotations_file="", min_pident=90, min_qcov=85):

    df = pd.read_csv(blast_outputs_file, delimiter="\t")
    if annotations_file != "":
        annotations = pd.read_csv(annotations_file)
        annotations[["cluster", "subcluster"]] = annotations["cluster"].str.split(" ", expand=True)
    df[["GCA_ID", "qseqid"]] = df["qseqid"].astype(str).str.split("|", expand=True) ## WGS sequences
    df[["gene", "sseqid"]] = df["sseqid"].str.split("|", expand=True)[[0,1]] ## reference sequences

    print("Number of total BLAST hits: " + str(len(df)) + "\n")
    print("Genes identified before filtering: " + str(len(df["gene"].unique())))
    min_pident = 90
    min_qcov = 85

    df["flo_cov"] = df["nident"] / df["length"] * 100 ## returns the coverage percentage of the portion of the query that was a match, as per Florentin's implementation in file 04

    df = df[(df["pident"] >= min_pident) & (df["flo_cov"] >= min_qcov)]
    print("Number of BLAST hits with % identity >= " + str(min_pident) + " and coverage >= " + str(min_qcov) + " : " + str(len(df)) + "\n")
    print("Genes identified after filtering: " + str(len(df["gene"].unique())))
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

    ## appending annotations / clusters for coloring
    if annotations_file != "":
        out_df = out_df.set_index("gene").join(annotations.set_index("gene"), lsuffix="_left", rsuffix="")

    return out_df