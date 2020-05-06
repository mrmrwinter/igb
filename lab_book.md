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

#-------------------------------

This is the point where this lab notebook was started.
The workflow runs at this point without error.
It produces:
1, Using code i have written - two pdfs of histograms of each sample. X axis is percent identity and Y axis is Frequency. One histo is just the second hit blast results, the other is the second hit blast results with all of the 100% hits removed.
2, Using Amirs Intragenomic blast script - two pngs of histograms of each sample. X axis is percent identity, im assuming Y is gene pairs or frequency? Again, one is the raw cdss, the other is the raw cdss with seqs that second hit at 100% removed.
