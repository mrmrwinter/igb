import intra_specific_blast as isb
from Bio import *
import io


cdss = snakemake.input[0]
cds_prefix = snakemake.params[0]

cds_prefix_isb = isb.IntraSpecificBlastCommandLine(cdss)
cds_prefix_isb.execute()
cds_prefix_percents = cds_prefix_isb.calculate_percents()

import matplotlib
import matplotlib.pyplot
from matplotlib import pyplot as plt
#plt.switch_backend('agg')
cds_prefix_hist = plt.hist(cds_prefix_percents, bins=200)
plt.plot(cds_prefix_hist[1][1:-1],cds_prefix_hist[0][:-1], color='red', linewidth=1)
plt.savefig(snakemake.output[0], dpi=900, frameon=False)
