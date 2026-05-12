import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

## importing the file
df = pd.read_csv("outputs/blast_outputs.tsv", delimiter="\t")
min_qcov = 30 ## filtering out the completely trash data
qcov_thresh = 75
pid_thresh = 85
df_filtered = df[df["qcovs"] > min_qcov]
x_axis = np.arange(0,100,0.1)
y_axis = np.arange(0,100,.1)


print(df_filtered["qcovs"].head())

sns.jointplot(data = df_filtered, x="qcovs", y="pident", kind="kde")
plt.xlabel("Query Coverage")
plt.xlim(min_qcov, 100)
plt.ylim(0, 100)
plt.vlines(qcov_thresh, 0, 100, linestyles="dashed")
plt.hlines(pid_thresh, 0, 100, linestyles="dashed")
plt.fill_between(x_axis, 100, pid_thresh, where = y_axis > qcov_thresh, alpha = 0.3, color = "green")
plt.ylabel("Percent Identity")
plt.title("BLAST Results", y=.9)
plt.show()