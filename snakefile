# ------------------------------------------
# MIKES INTRAGENOMIC CDS BLASTING SNAKE
# ------------------------------------------

# have fasta with 1000s of cds's in
SAMPLE="data/input/{sample}.cds_nt.fa"

rule all:
    input:
        "outputs/plot.png"
# make database of cds fasta
rule make_blast_database:
    input:
        "data/database/{sample}.cds_nt.fa"
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
        db = "data/input/{sample}_cds",
        query = "data/database/{sample}.cds_nt.fa"
    output:
        "outputs/blast_results"
    params:
        "data/database/{sample}_cds"
    threads:
        8
    shell:
        "blastn -q {input} -db {params} -o {output}"


# -------------------------------------------
import pandas as pd

def read_blast_output(output):
    """Reads BLAST output (outfmt 6) and returns a pandas dataframe."""

	return pd.read_csv(output,
                       sep="\t",
                       names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"],
                       index_col="qseqid")
# -----------------------------------------------



# take second top hit
rule take_second_hit:
# needs to parse blast output and take the second top hit row and write it to a table
    input:
        "outputs/blast_results"
    output:
        "outputs/second_hits"
    run:
        """
        for blast_result in blast_results
            do
                # take row with second top hit
                #append it to results_table.tsv
        """


# take percent ident of all seconds hits
rule take_percent_ident:
# this will take the column with all of the percent identities and either add them to a df or a table
    input:
        "outputs/second_hits"
    output:
        "outputs/percent_idents"
    shell:
        # sed/awk row with idents in and print to a file


#plot
rule plot:
    input:
        "outputs/percent_idents"
    output:
        "outputs/plots/{sample}.png"
    shell:
        "Rscript plot"
# this will plot the percent identities as a histogram
