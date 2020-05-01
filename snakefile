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
        expand("outputs/{sample}_blastResults", sample=SAMPLE),
        expand("outputs/{sample}_secondHits", sample=SAMPLE),
        expand("outputs/{sample}_percentIdents", sample=SAMPLE),
        expand("outputs/{sample}_percentIdents", sample=SAMPLE),
        expand("outputs/{sample}/{sample}.pdf", sample=SAMPLE),
        expand("outputs/{sample}_jupInNo100", sample=SAMPLE)

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
        "outputs/{sample}_blastResults"
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
rule remove_hundreds:
    input:
        "outputs/{sample}_blastResults"
    output:
        "outputs/{sample}_blastResultsNoHundred"
    shell:
        "grep -v '100.00' {input} > {output}"   # find a way to do this in python
                                                # the blast results should be moved to a pandas frame before this step
#
# # take second top hit
rule take_second_hit:
# needs to parse blast output and take the second top hit row and write it to a table
    input:
        "outputs/{sample}_blastResultsNoHundred"
    output:
        "outputs/{sample}_secondHits"
    shell:
      "awk 'NR % 2 == 0' {input} > {output}"
#




rule make_input_for_igbpynb:
# needs to parse blast output and take the second top hit row and write it to a table
    input:
        "outputs/{sample}_secondHits"
    output:
        "outputs/{sample}_cdssForJupyter"
    shell:
      "grep -v '100.00' {input} | cut -f1 > {output}"
#
rule make_input_for_igbpynb_2:
    input:
        headers = "outputs/{sample}_cdssForJupyter",
        cdss = "data/input/{sample}.cds_nt.fa"
    output:
        "outputs/{sample}_jupInNo100"
    shell:
        "seqtk subseq {input.cdss} {input.headers} > {output}"
#
# # take percent ident of all seconds hits
rule take_percent_ident:
# this will take the column with all of the percent identities and either add them to a df or a table
    input:
        "outputs/{sample}_secondHits"
    output:
        "outputs/{sample}_percentIdents"
    shell:
        "awk '{{print $3}}' {input} | grep -v '100.000' > {output}"


# #plot
rule plot:
    input:
        "outputs/{sample}_percentIdents"
    output:
        "outputs/{sample}/{sample}.pdf"
    script:
        "scripts/pdplot.py"
