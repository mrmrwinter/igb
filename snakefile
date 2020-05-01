# ------------------------------------------
# MIKES INTRAGENOMIC CDS BLASTING SNAKE
# ------------------------------------------

import os

# input data
## samples = {f[:-8] for f in os.listdir("data/input") if f.endswith(".cds_nt.fa")}

# have fasta with 1000s of cds's in
SAMPLE, = glob_wildcards("data/input/{sample}.cds_nt.fa")

rule all:
    input:
        expand("data/database/{sample}_cds.nhr", sample=SAMPLE),
        expand("outputs/{sample}_blastResults", sample=SAMPLE),
        expand("outputs/{sample}_secondHits", sample=SAMPLE),
        expand("outputs/{sample}_percentIdents", sample=SAMPLE),
        expand("outputs/{sample}_percentIdents", sample=SAMPLE),
        expand("outputs/{sample}.pdf", sample=SAMPLE)

# make database of cds fasta
rule make_blast_database:
    input:
        "data/input/{sample}.cds_nt.fa"
    output:
        nhr = "data/database/{sample}_cds.nhr",
        nin = "data/database/{sample}_cds.nin",
        nsq = "data/database/{sample}_cds.nsq"
    params:
        "data/database/{sample}_cds"
    priority:
        50
    threads:
        8
    shell:
        "makeblastdb -in {input} -out {params} -dbtype nucl"

# Blast the cds msa against blast_database
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
#
# # take second top hit
rule take_second_hit:
# needs to parse blast output and take the second top hit row and write it to a table
    input:
        "outputs/{sample}_blastResults"
    output:
        "outputs/{sample}_secondHits"
    shell:
      "awk 'NR % 2 == 0' {input} > {output}"
#
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
        "outputs/{sample}.pdf"
    script:
        "scripts/pdplot.py"
