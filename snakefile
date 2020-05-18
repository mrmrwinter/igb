# ------------------------------------------
# MIKES INTRAGENOMIC CDS BLASTING SNAKE
# ------------------------------------------

# Place input files in data/input/
# Ensure they end in .cds_nt.fa
# Ensure headers of fastas contain no spaces

import os  # imports os

# input data
## samples = {f[:-8] for f in os.listdir("data/input") if f.endswith(".cds_nt.fa")}  # this is a potential way to do iteration

SAMPLE, = glob_wildcards("data/input/{sample}.cds_nt.fa")  # globs the sample names of all files

rule all_outputs:  # target rule
    input:
        expand("data/database/{sample}_cds.nhr", sample=SAMPLE),  # ensures creation of database
        expand("outputs/{sample}/blastResults", sample=SAMPLE),  # results from blast search
        expand("outputs/{sample}/secondHits", sample=SAMPLE),  # second hits from blast
        expand("outputs/{sample}/percentIdents", sample=SAMPLE),  # percent idents
        expand("outputs/{sample}/percentIdentsNo100", sample=SAMPLE),  # percent idents with 100s removed
        expand("results/{sample}/{sample}.pdf", sample=SAMPLE),  # snakemake plot
        expand("results/{sample}/{sample}No100.pdf", sample=SAMPLE),  # snakemake plot with 100s removed
        expand("outputs/{sample}/jupInNo100.cds_nt.fa", sample=SAMPLE),  # file to go into the notebook
        #expand("results/{sample}/pyplot.png", sample=SAMPLE),
        expand("results/{sample}/igbpyOut.png", sample=SAMPLE),  # notebook output
        expand("results/{sample}/igbpyOutNo100.png", sample=SAMPLE)  # notebook output from feeding percents in halfway through

# -----------------------------------------
# DATABASE CREATION

rule make_blast_database:  # Rule to make database of cds fasta
    input:
        "data/input/{sample}.cds_nt.fa" # input to the rule
    output:
        nhr = "data/database/{sample}_cds.nhr",   # all outputs expected from the rule
        nin = "data/database/{sample}_cds.nin",
        nsq = "data/database/{sample}_cds.nsq"
    params:
        "data/database/{sample}_cds"  # prefix for the outputs, required by the command
    priority:
        50   # priority of rule. makes rule be ran first regardless of I/O
    threads:
        8  # number of threads to use
    shell:  # shell command for the rule
        "makeblastdb \
        -in {input} \
        -out {params} \
        -dbtype nucl"  # the database type


# -------------------------------------------------------------
# BLASTING

rule blastn:  # Blast the cds multiFASTA against blast_database
    input:
        db = "data/database/{sample}_cds.nhr",  # blast database
        query = "data/input/{sample}.cds_nt.fa"  # input query
    output:
        "outputs/{sample}/blastResults"  # outputs a table of blast results
    params:
        "data/database/{sample}_cds"  # prefix for the command
    threads:
        8
    shell:
        "blastn \
        -db {params} \
        -query {input.query} \
        -max_target_seqs 2 \
        -outfmt '6 qseqid sseqid pident evalue' \
        -out {output}"
        # use blastn
        #  which database to use
        #   the query input
        #    the amount of hits to return
        #     the output format - .tsv with query seqid, hit seqid, percent identity, and evalue
        #      the output destination


# ---------------------------------------------------------
# DATA TRANSFORMATION

rule take_second_hit:  # parses blast output, takes the second top hit row and writes it to a table
    input:
        "outputs/{sample}/blastResults"  # input is the blast results
    output:
        "outputs/{sample}/secondHits"  # outputs a .tsv of only the second top hits
    shell:
      "awk 'NR%2==0' {input} > {output}"  # takies second row from a table and writes to a new table


rule remove_hundreds:  # removes any results in secondHits that hit at 100%
    input:
        "outputs/{sample}/secondHits"
    output:
        "outputs/{sample}/secondHitsNo100"
    shell:
        "grep -v '100.00' {input} > {output}"   # do this in python
                                                # the blast results should be moved to a pandas frame before this step


rule take_percent_ident:  # take percent ident of all seconds hits and either add them to a new table
    input:
        "outputs/{sample}/secondHits"  # secondHit table
    output:
        "outputs/{sample}/percentIdents"  # table of only percent identities
    shell:
        "awk '{{print $3}}' {input} > {output}"  # takes the percent identity column and writes to new table


rule take_percent_ident_no100:  # this does the same as the rule above but does not copy over hits of 100%
    input:
        "outputs/{sample}/secondHitsNo100"
    output:
        "outputs/{sample}/percentIdentsNo100"
    shell:
        "awk '{{print $3}}' {input} > {output}"  # takes percent ident column and writes to new table


# -------------------------------------------------------------------
# PLOTTING

rule plot_no100:  # plots the percent identities with 100%s removed as a histogram
    input:
        "outputs/{sample}/percentIdentsNo100"  # percent ident table with 100% hits removed
    output:
        "results/{sample}/{sample}No100.pdf"  # outputs a histogram in pdf
    script:
        "scripts/pdplot.py"  # runs from a python script


rule plot:  # plots the percent identities as a histogram
    input:
        "outputs/{sample}/percentIdents"  # percentIdent table
    output:
        "results/{sample}/{sample}.pdf"  # outputs a histogram in pdf format
    script:
        "scripts/pdplot.py"  # runs a python script that inputs to pandas and creates a plot


rule igb:  # runs a script written from Amirs notebooks
    input:
        "data/input/{sample}.cds_nt.fa"  # input is the "raw" cds data
    output:
        "results/{sample}/igbpyOut.png"  # outputs a smoothed histogram in png
    params:
        "{sample}"  # prefix for the script
    script:
        "igbsm.py"  # runs from a python script

rule igbNo100:  # runs a script written from Amirs notebooks
    input:
        "data/input/{sample}.cds_nt.fa"  # input is the "raw" cds data
    output:
        "results/{sample}/igbpyOut.png"  # outputs a smoothed histogram in png
    params:
        "{sample}/"  # prefix for the script
    script:
        "igbsm.py"  # runs from a python script


### code below for overlayed line hist

# plt.plot(MfloJB5_hist[1][1:-1],MfloJB5_hist[0][:-1], color='red', linewidth=2)
# plt.plot(mhap_63.nh_hist[1][1:-1],mhap_63.nh_hist[0][:-1], color='grey', linewidth=2)
# plt.plot(MincW1_hist[1][1:-1],MincW1_hist[0][:-1], color='green', linewidth=2)
# plt.plot(MjavVW4_hist[1][1:-1],MjavVW4_hist[0][:-1], color='green', linewidth=2)
# plt.savefig('hap_inc_flo_blast.png', dpi=900, frameon=False)


# for worm in {input}:


# -------------------------------------------
# import pandas as pd
#
# def read_blast_output(output):
#     #"""Reads BLAST output (outfmt 6) and returns a pandas dataframe."""
#
# 	return pd.read_csv(output,
#                        sep="\t",
#                        names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"],
#                        index_col="qseqid")
# -----------------------------------------------

#

## this was old code for making an inpit for the notebook
## fairly sure this is redundant now
# rule make_input_for_igbpy:
# # needs to parse blast output and take the second top hit row and write it to a table
#     input:
#         "outputs/{sample}/secondHits"
#     output:
#         "outputs/{sample}/cdssForJupyter"
#     shell:
#       "grep -v '100.00' {input} | cut -f2 > {output}"
# #
# rule make_input_for_igbpy_2:
#     input:
#         headers = "outputs/{sample}/cdssForJupyter",
#         cdss = "data/input/{sample}.cds_nt.fa"
#     output:
#         "outputs/{sample}/jupInNo100.cds_nt.fa"
#     shell:
#         "seqtk subseq {input.cdss} {input.headers} > {output}"
# rule igb_No100:  # runs the plotting section of Amirs notebooks, inputting percentIdents with 100%s removed
#     input:
#         "outputs/{sample}/jupInNo100.cds_nt.fa"
#     output:
#         "results/{sample}/igbpyOutNo100.png"
#     params:
#         "{sample}"
#     script:
#         "igbsm.py"
