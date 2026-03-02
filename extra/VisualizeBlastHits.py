import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

## importing the file
df = pd.read_csv("outputs/blast_outputs.tsv", delimiter="\t")


sns.jointplot(data = df, x="qcovs", y="pident", kind="kde")
plt.xlabel("Query Coverage")
plt.ylabel("Percent Identity")
plt.title("BLAST Results", y=.9)
plt.show()