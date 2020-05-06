# Lab book for Intragenomic blast

## Experiments and progress of the intragenomic blast workflow will be recorded here

### Raw cds data are in data/raw/
### Inputs for the workflow are in data/input/
### Outputted plots are in results/{sample}
### Intermediary files are in outputs/{sample}
### Remote repository is https://github.com/mrmrwinter/igb

### The workflow is ran using
`snakemake --cores 8 -s snakefile`

### The lab book will be in reverse chronological order from here

## 6/5/20

Cleaned up directories.

Started adding the rest of Amirs IGB notebook to the snakemake workflow.
Amirs notebook can be found here: https://tinyurl.com/y8k62mw6

Added a cell to igb.ipynb with the code:

`plt.plot(minc_hist[1][1:-1],minc_hist[0][:-1], color='red', linewidth=2)
plt.plot(hapnh_hist[1][1:-1],hapnh_hist[0][:-1], color='grey', linewidth=2)
plt.plot(mflo_hist[1][1:-1],mflo_hist[0][:-1], color='green', linewidth=2)
plt.plot(jav_hist[1][1:-1],jav_hist[0][:-1], color='blue', linewidth=2)
plt.savefig('hap_inc_flo_blast.png', dpi=900, frameon=False)`

this plots the second hit percent idents of incognita, floridensis, and haplanaria.
Will need the cdss of the uncollapsed haplanaria assembly as well as cdss from hapla for comparison.
Javanica (MjavVW4) not working in this script for some reason. Returns error: `File is not accessible:  "data/input/MjavVW4.cds_nt.fa.98\"`
This was due to the path being defined incorrectly in an earlier cell.

Moved onto plotting the smoothed seaborn histogram with a legend.

Copied the code from Amirs notebook: https://tinyurl.com/y8k62mw6
Installed necessary dependencies into IGB environment:
`conda activate IGB`
`conda install seaborn`

Reproduced Amirs plot with the following script:

`fig, ax = plt.subplots(figsize=(10,8))
bins = []

colors = {'jav':'blue','mflo':'black','hapnh':'green',
          'minc':'red'}

# sns.distplot(a,
#             hist=False,
#             color="blue",
#             kde_kws={'shade': False, 'lw':w},
#             #rug=True,
#             )

sns.distplot(b,
            hist=False,
            color="blue",
            kde_kws={'shade': False, 'lw':w},
            #rug=True,
            )

sns.distplot(c,
            hist=False,
            color="red",
            kde_kws={'shade': False, 'lw':w},
            #rug=True,
            )

sns.distplot(d,
            hist=False,
            color="green",
            kde_kws={'shade': False, 'lw':w},
            #rug=True,
            )

ax.set_xlim(78,100)

plt.xlabel('Percent Identity', fontsize=18)
plt.ylabel('Prop. Gene Pairs', fontsize=18)

# Plot the legend
Mhaplanaria_line = mlines.Line2D([0,1], [0,0], color='green',
                          linewidth=3, label='M. haplanaria')
Minc_line = mlines.Line2D([0,1], [0,0], color='blue',
                          linewidth=3, label='M. incognita W1')
Mflo_line = mlines.Line2D([0,1], [0,0], color='red',
                          linewidth=3, label='M. floridensis SJF1/JB5')

plt.legend(handles=[Minc_line,Mflo_line,Mhaplanaria_line],
           loc=2, fontsize = 16)`

Javanica wont work in this cell/script. Gives error: `TypeError: len() of unsized object`


### To do next
Get the latest python plotting scripts working in the snakemake workflow.


#-------------------------------

This is the point where this lab notebook was started.
The workflow runs at this point without error.
It produces:
1, Using code i have written - two pdfs of histograms of each sample. X axis is percent identity and Y axis is Frequency. One histo is just the second hit blast results, the other is the second hit blast results with all of the 100% hits removed.
2, Using Amirs Intragenomic blast script - two pngs of histograms of each sample. X axis is percent identity, im assuming Y is gene pairs or frequency? Again, one is the raw cdss, the other is the raw cdss with seqs that second hit at 100% removed.
