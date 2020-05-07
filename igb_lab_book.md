# Lab book for Intragenomic blast

### Experiments and progress of the intragenomic blast workflow will be recorded here

 Raw cds data are in data/raw/

 Inputs for the workflow are in data/input/

 Outputted plots are in results/{sample}

 Intermediary files are in outputs/{sample}

 Remote repository is https://github.com/mrmrwinter/igb

 The workflow is ran using

`snakemake --cores 8 -s snakefile`

#### The lab book will be in reverse chronological order from here




### 7/5/20

#### Things from before
Get the latest python plotting scripts working in the snakemake workflow.

Get it ran with uncollapsed assembly.

Sort out the removal of the 100s




#### things done

Sorted removal of 100% hits in M. haplanaria by inputting the `outputs/mhap_63_nh/percentIdents` to amirs notebook after the blast stage

Changed paths in Amirs notebook to fit the data structure

Can now replicate Amirs igb plot exactly using scripts in the `intraspecific_blast.ipynb`




### 6/5/20

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

Reproduced Amirs plot

Javanica wont work in this cell/script. Gives error: `TypeError: len() of unsized object`

Also hap has a large spike on 100. this could be due to the assembly being collapsed.

Attempts to remove cdss that have a second hit of 100 have failed so far: see 'rule make_input_for_igbpy' in the snakefile


### To do next
Get the latest python plotting scripts working in the snakemake workflow.

Get it ran with uncollapsed assembly.

Sort out the removal of the 100s


#-------------------------------

This is the point where this lab notebook was started.

The workflow runs at this point without error.

It produces:

1, Using code i have written - two pdfs of histograms of each sample. X axis is percent identity and Y axis is Frequency. One histo is just the second hit blast results, the other is the second hit blast results with all of the 100% hits removed.

2, Using Amirs Intragenomic blast script - two pngs of histograms of each sample. X axis is percent identity, im assuming Y is gene pairs or frequency? Again, one is the raw cdss, the other is the raw cdss with seqs that second hit at 100% removed.
