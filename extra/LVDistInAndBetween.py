## sanity check test that takes each gene and calculates the Levenstein distance 
## within a given variant and within all other variants to ensure that variant labels are sensible
## This produces a modified dataframe that 
import pandas as pd
from DFfromBLAST import create_df_from_blast
import Levenshtein as lv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

blast_outputs_file = "outputs/test.tsv"
annotations_file = "extra/gene_annotations.csv"
min_pident = 90
min_qcov = 85

def create_LV_dist_df(blast_outputs_file, annotations_file = "", min_pident = 90, min_qcov=85, create_plots = False):
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

    ## filtering to just BLP
    if create_plots == True:
        out_df = out_df[out_df["cluster"] == "blp"]
        
        fig, ax = plt.subplots()

        print(out_df.function)
        sns.scatterplot(data=out_df, x="ingroup_mean_ratio", y="outgroup_mean_ratio", hue="function")#, style = "function")
        ax.set_xlabel("mean ingroup similarity ratio")
        ax.set_ylabel("mean outgroup similarity ratio")
        scatter_lim_x = ax.get_xlim()
        scatter_lim_y = ax.get_ylim()
        plt.plot([0,1], [0,1], linestyle="dashed")
        ax.set_xlim(scatter_lim_x)
        ax.set_ylim(scatter_lim_y)
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1,1))
        plt.show()

    return out_df