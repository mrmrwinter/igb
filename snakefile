# ------------------------------------------
# MIKES INTRAGENOMIC CDS BLASTING SNAKE
# ------------------------------------------


import os   #imports os

# input data
## samples = {f[:-8] for f in os.listdir("data/input") if f.endswith(".cds_nt.fa")}
# this is a potential way to do iteration


# have fasta with 1000s of cds's in
SAMPLE, = glob_wildcards("data/input/{sample}.cds_nt.fa")
# globs the sample names of all files


rule all_outputs:
    input:
        expand("data/database/{sample}_cds.nhr", sample=SAMPLE),
        expand("outputs/{sample}/blastResults", sample=SAMPLE),
        expand("outputs/{sample}/secondHits", sample=SAMPLE),
        expand("outputs/{sample}/percentIdents", sample=SAMPLE),
        expand("outputs/{sample}/percentIdents", sample=SAMPLE),
        expand("results/{sample}/{sample}.pdf", sample=SAMPLE),
        expand("outputs/{sample}/jupInNo100.cds_nt.fa", sample=SAMPLE),
        #expand("results/{sample}/pyplot.png", sample=SAMPLE),
        expand("results/{sample}/igbpyOut.png", sample=SAMPLE),
        expand("results/{sample}/igbpyOutNo100.png", sample=SAMPLE)

# make database of cds fasta
rule make_blast_database:    # name of thje rule
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
    shell:
        "makeblastdb -in {input} -out {params} -dbtype nucl"  # the database type


# Blast the cds multiFASTA against blast_database
rule blastn:
    input:
        db = "data/database/{sample}_cds.nhr",
        query = "data/input/{sample}.cds_nt.fa"
    output:
        "outputs/{sample}/blastResults"
    params:
        "data/database/{sample}_cds"
    threads:
        8
    shell:
        "blastn -db {params} -query {input.query} -max_target_seqs 2 -outfmt '6 qseqid sseqid pident evalue' -out {output}"


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

#
# # take second top hit
rule take_second_hit:
# needs to parse blast output and take the second top hit row and write it to a table
    input:
        "outputs/{sample}/blastResults"
    output:
        "outputs/{sample}/secondHits"
    shell:
      "awk 'NR % 1 == 0' {input} > {output}"
#

rule remove_hundreds:
    input:
        "outputs/{sample}/secondHits"
    output:
        "outputs/{sample}/secondHitsNoHundred"
    shell:
        "grep -v '100.00' {input} > {output}"   # find a way to do this in python
                                                # the blast results should be moved to a pandas frame before this step


rule make_input_for_igbpy:
# needs to parse blast output and take the second top hit row and write it to a table
    input:
        "outputs/{sample}/secondHits"
    output:
        "outputs/{sample}/cdssForJupyter"
    shell:
      "grep -v '100.00' {input} | cut -f2 > {output}"
#
rule make_input_for_igbpy_2:
    input:
        headers = "outputs/{sample}/cdssForJupyter",
        cdss = "data/input/{sample}.cds_nt.fa"
    output:
        "outputs/{sample}/jupInNo100.cds_nt.fa"
    shell:
        "seqtk subseq {input.cdss} {input.headers} > {output}"
#
# # take percent ident of all seconds hits
rule take_percent_ident:
# this will take the column with all of the percent identities and either add them to a df or a table
    input:"outputs/{sample}/secondHits"
    output:
        "outputs/{sample}/percentIdents"
    shell:
        "awk '{{print $3}}' {input} | grep -v '100.000' > {output}"


# #plot
rule plot:
    input:
        "outputs/{sample}/percentIdents"
    output:
        "results/{sample}/{sample}.pdf"
    script:
        "scripts/pdplot.py"


# #### python method
#
# rule plotting:
#     input:
#         "outputs/{sample}/percentIdents",
#     output:
#         "results/{sample}/pyplot.png"
#     params:
#         "{sample}"
#     run:
#         """
#         import matplotlib
#         {params}_hist = plt.hist({input}, bins=200)
#         plt.plot({params}_hist[1][1:-1],{params}_hist[0][:-1], color='red', linewidth=2)
#         plt.savefig('{output}', dpi=900, frameon=False)
#         """

rule igb:
    input:
        "data/input/{sample}.cds_nt.fa"
    output:
        "results/{sample}/igbpyOut.png"
    params:
        "{sample}"
    script:
        "igbsm.py"

rule igb_No100:
    input:
        "outputs/{sample}/jupInNo100.cds_nt.fa"
    output:
        "results/{sample}/igbpyOutNo100.png"
    params:
        "{sample}"
    script:
        "igbsm.py"
