from Bio import *
import io


idents = snakemake.input[0]
idents_prefix = snakemake.params[0]

import matplotlib
import matplotlib.pyplot
from matplotlib import pyplot as plt
#plt.switch_backend('agg')
idents_prefix_hist = plt.hist(idents_prefix_percents, bins=200)
plt.plot(idents_prefix_hist[1][1:-1],idents_prefix_hist[0][:-1], color='red', linewidth=1)
plt.savefig(snakemake.output[0], dpi=900, frameon=False)
