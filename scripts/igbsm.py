import intra_specific_blast2 as isb
import Bio
import io


cdss = snakemake.input[0]
cds_prefix = snakemake.params[0]

cds_prexix_isb = isb.IntraSpecificBlastCommandLine(cdss)
cds_prexix_isb.execute()
cds_prexix_percents = cds_prexix_isb.calculate_percents()

from reprophylo import *
import matplotlib.pyplot
plt.switch_backend('agg')
cds_prexix_hist = plt.hist(cds_prexix_percents, bins=200)
plt.plot(cds_prexix_hist[1][1:-1],cds_prexix_hist[0][:-1], color='red', linewidth=2)
plt.savefig(snakemake.output[0], dpi=900, frameon=False)
