################################################################################
# Project: Training project
# Part: Composition analysis
# Step: Taxonomic profiling with Kraken2
#
# Dependent on:
#   - PREP_remove_hostDNA.Snakefile
#
# Alex Huebner, 18/04/23
################################################################################

import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("03-data/processed_data/{sample}_1.fastq.gz")
################################################################################

KMERLENS = [50, 75, 100, 150, 200, 250, 300]

rule all:
    input:
        "05-results/COMP_Kraken2_standardDB.tsv.gz",
        "05-results/COMP_Bracken_Kraken2_standardDB.tsv.gz"

#### Kraken2 ###################################################################

rule kraken2_download_db:
    output:
        "03-data/refdbs/kraken2_standard_20221209/hash.k2d"
    message: "Download and extract Kraken2 database"
    params:
        url = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20221209.tar.gz",
        tarball = "tmp/k2_standard_20221209.tar.gz",
        basedir = "03-data/refdbs"
    shell:
        """
        wget -O {params.tarball} {params.url}
        mkdir -p {params.basedir}/$(basename {params.tarball} .tar.gz)
        tar xvf {params.tarball} -C {params.basedir}/$(basename {params.tarball} .tar.gz)
        """

rule kraken2:
    input:
        "03-data/refdbs/kraken2_standard_20221209/hash.k2d"
    output:
        "04-analysis/kraken2/{sample}.report.txt"
    message: "Screen against Kraken2 standard database: {wildcards.sample}"
    conda: "ENVS_Kraken2_Bracken.yaml"
    resources:
        mem = 80,
        cores = 16
    params:
        db = "03-data/refdbs/kraken2_standard_20221209",
        pe1 = "03-data/processed_data/{sample}_1.fastq.gz",
        pe2 = "03-data/processed_data/{sample}_2.fastq.gz"
    threads: 16
    shell:
        """
         kraken2 --db {params.db} \
             --threads {threads} \
             --output - \
             --report {output} \
             --confidence 0.15 \
             --paired --use-names --gzip-compressed \
             {params.pe1} {params.pe2}
        """

rule summarise_kraken2:
    input:
        expand("04-analysis/kraken2/{sample}.report.txt", sample=SAMPLES)
    output:
        temp("05-results/COMP_Kraken2_standardDB.tsv")
    message: "Summarise the Kraken2 profiles using Taxpasta"
    conda: "ENVS_taxpasta.yaml"
    resources:
        mem = 8
    shell:
        """
        taxpasta merge -p kraken2 -o {output} --long {input}
        """

rule reformat_taxpasta:
    input:
        "05-results/COMP_Kraken2_standardDB.tsv"
    output:
        "05-results/COMP_Kraken2_standardDB.tsv.gz"
    message: "Re-format Taxpasta output"
    resources:
        mem = 8
    run:
        res = pd.read_csv(input[0], sep="\t")
        res['sample'] = res['sample'].str.replace(".report", "")
        res[['sample', 'taxonomy_id', 'count']] \
            .to_csv(output[0], sep="\t", index=False, compression="gzip")

################################################################################

#### Bracken ###################################################################

rule bracken:
    input:
        "04-analysis/kraken2/{sample}.report.txt"
    output:
        txt = "04-analysis/bracken/{sample}.bracken_{l}bp.txt"
    message: "Infer relative abundances with Bracken: {wildcards.sample} assuming read length {wildcards.l} bp"
    conda: "ENVS_Kraken2_Bracken.yaml"
    resources:
        mem = 2,
        cores = 1
    params:
        db = "03-data/refdbs/kraken2_standard_20221209"
    benchmark: "tmp/bracken.{sample}_{l}.benchmark.txt"
    threads: 1
    shell:
        """
        bracken -d {params.db} \
                -i {input} \
                -o {output.txt} \
                -r {wildcards.l} \
                -l S
        """

rule summarise_bracken:
    input:
        expand("04-analysis/bracken/{sample}.bracken_{l}bp.txt", sample=SAMPLES, l=KMERLENS)
    output:
        "05-results/COMP_Bracken_Kraken2_standardDB.tsv.gz"
    message: "Summarise the taxonomic profiles produced by Bracken"
    resources:
        mem = 8,
        cores = 1
    run:
        reports = pd.concat([pd.read_csv(fn, sep="\t") \
                             .assign(filename=os.path.basename(fn))
                             for fn in input])
        reports['sample'] = reports['filename'].str.split(".").str[0]
        reports['kmerlength'] = reports['filename'].str.extract(r'[ES]RS[0-9]+.bracken_([0-9]+)bp.txt').astype(int)
        reports[['sample', 'kmerlength', 'name', 'taxonomy_id', 'new_est_reads', 'fraction_total_reads']] \
            .rename({'name': 'taxon',
                     'new_est_reads': 'count',
                     'fraction_total_reads': 'relAb'}, axis=1) \
            .to_csv(output[0], sep="\t", index=False, compression="gzip")

################################################################################
