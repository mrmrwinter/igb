import intra_specific_blast2 as isb
import Bio
import io


cdss = {input}
{params}_isb = isb.IntraSpecificBlastCommandLine(cdss)
{params}_isb.execute()
{params}_percents = {params}_isb.calculate_percents()

from reprophylo import *
%matplotlib inline
plt.switch_backend('agg')
{params}_hist = plt.hist({params}_percents, bins=200)
plt.plot({params}_hist[1][1:-1],{params}_hist[0][:-1], color='red', linewidth=2)
plt.savefig('{output}', dpi=900, frameon=False)
