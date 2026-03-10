## sanity check test that takes each gene and calculates the Levenstein distance 
## within a given variant and within all other variants to ensure that variant labels are sensible
## This produces a modified dataframe that 
import pandas as pd
from DFfromBLAST import create_df_from_blast
import Levenshtein as lv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

## need a way to access my helper functions from other folder
import sys
import os
cwd = os.getcwd()
new_path = cwd.split("/")
new_path.append("scripts") 
sys.path.append("/".join(new_path))
from PythonHelperFuncs import create_LV_dist_df


blast_outputs_file = "outputs/blast_out.tsv"
annotations_file = "extra/gene_annotations.csv"
min_pident = 90
min_qcov = 85
plot_target_cluster = "blp"

out_df = create_LV_dist_df(blast_outputs_file, annotations_file, min_pident, min_qcov, verbose=True)

if plot_target_cluster != "":
    out_df = out_df[out_df["cluster"] == plot_target_cluster]
fig, ax = plt.subplots()
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
plt.savefig("blp_cluster_sequence_comparisons.png")