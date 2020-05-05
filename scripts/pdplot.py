import pandas as pd
import matplotlib

df = pd.read_csv(snakemake.input[0], sep='\t')
ax = df.plot.hist(bins=100, alpha=0.5)
fig = ax.get_figure()
fig.savefig(snakemake.output[0])
