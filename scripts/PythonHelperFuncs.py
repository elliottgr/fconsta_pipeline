import pandas as pd
import numpy as np
import Levenshtein as lv

def create_df_from_blast(blast_outputs_file, annotations_file="", min_pident=90, min_qcov=85, verbose = False):

    df = pd.read_csv(blast_outputs_file, delimiter="\t")
    if annotations_file != "":
        annotations = pd.read_csv(annotations_file)
        annotations[["cluster", "subcluster"]] = annotations["cluster"].str.split(" ", expand=True)
    df[["GCA_ID", "qseqid"]] = df["qseqid"].astype(str).str.split("|", expand=True) ## WGS sequences
    df[["gene", "sseqid"]] = df["sseqid"].str.split("|", expand=True)[[0,1]] ## reference sequences

    if verbose == True:
        print("Number of total BLAST hits: " + str(len(df)) + "\n")
        print("Genes identified before filtering: " + str(len(df["gene"].unique())))

    df["flo_cov"] = df["nident"] / df["length"] * 100 ## returns the coverage percentage of the portion of the query that was a match, as per Florentin's implementation in file 04

    df = df[(df["pident"] >= min_pident) & (df["flo_cov"] >= min_qcov)]
    if verbose == True:
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

def create_LV_dist_df(blast_outputs_file, annotations_file = "", min_pident = 90, min_qcov=85, create_plots = False, plot_target_cluster = "", verbose = False):
    df = create_df_from_blast(blast_outputs_file, annotations_file, min_pident, min_qcov)

    ## variant with only blp data
    blp_df = df[df["cluster"] == "blp"].groupby("gene")

    ## creating arrays for a tidy df later
    gene_id = []
    variant_id = []
    functions = []
    clusters = []

    ingroup_means = []
    ingroup_medians = []
    ingroup_means_ratio = []
    ingroup_medians_ratio = []
    outgroup_means = []
    outgroup_medians = []
    outgroup_means_ratio = []
    outgroup_medians_ratio = []

    ## iterate over each gene
    for gene, group in blp_df:
        ## per variant
        for variant in group["variant"].unique():
            if verbose == True:
                print("Now comparing " + str(gene) + " variant ID " + str(variant))
            within_group = group[group["variant"] == variant]
            out_group = group[group["variant"] != variant]

            ingroup_lvs = []
            outgroup_lvs = []
            ingroup_lvs_ratio = []
            outgroup_lvs_ratio = []
            for seq_foc in within_group.sequences:
                for seq_tar in within_group.sequences:
                    ingroup_lvs.append(lv.distance(seq_foc, seq_tar))
                    ingroup_lvs_ratio.append(lv.ratio(seq_foc, seq_tar))
                for seq_tar in out_group.sequences:
                    outgroup_lvs.append(lv.distance(seq_foc, seq_tar))
                    outgroup_lvs_ratio.append(lv.ratio(seq_foc, seq_tar))


            gene_id.append(gene)
            variant_id.append(variant)
            functions.append(group["functions"].unique()[0])
            clusters.append(group["cluster"].unique()[0])

            ingroup_means.append(np.mean(ingroup_lvs))
            ingroup_medians.append(np.median(ingroup_lvs))
            ingroup_means_ratio.append(np.mean(ingroup_lvs_ratio))
            ingroup_medians_ratio.append(np.median(ingroup_lvs_ratio))

            outgroup_means.append(np.mean(outgroup_lvs))
            outgroup_medians.append(np.median(outgroup_lvs))
            outgroup_means_ratio.append(np.mean(outgroup_lvs_ratio))
            outgroup_medians_ratio.append(np.median(outgroup_lvs_ratio))

    out_df = pd.DataFrame(data={"gene":gene_id, "variant":variant_id, "function":functions, "cluster":clusters,
                                "ingroup_mean_ratio":ingroup_means_ratio, "outgroup_mean_ratio":outgroup_means_ratio, 
                                "ingroup_median_ratio":ingroup_medians_ratio, "outgroup_median_ratio":outgroup_medians_ratio,
                                "ingroup_mean":ingroup_means, "outgroup_mean":outgroup_means, 
                                "ingroup_median":ingroup_medians, "outgroup_median":outgroup_medians})
    return out_df