################################################################################
# Project: Training project
# Part: Composition analysis
# Step: Taxonomic profiling with Kraken2 on Zape2
#
# Alex Huebner, 19/04/23
################################################################################

import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES, =  glob_wildcards("../week2/eager_Zape2/mapping/bwa/{sample}_PE.mapped.bam")
################################################################################

rule all:
    input:
        expand("04-analysis/kraken2/{sample}.report.txt", sample=SAMPLES)
    
rule bam2fq:
    output:
        pipe("tmp/kraken2/{sample}.fastq")
    message: "Extract unmapped reads and convert into FastQ: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    params:
        bam = "../week2/eager_Zape2/mapping/bwa/{sample}_PE.mapped.bam"
    shell:
        """
        samtools view -uh -f 4 {params.bam} | samtools fastq - > {output}
        """

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
        db = "03-data/refdbs/kraken2_standard_20221209/hash.k2d",
        fq = "tmp/kraken2/{sample}.fastq"
    output:
        "04-analysis/kraken2/{sample}.report.txt"
    message: "Screen against Kraken2 standard database: {wildcards.sample}"
    conda: "ENVS_Kraken2_Bracken.yaml"
    resources:
        mem = 100,
        cores = 16
    params:
        db = "03-data/refdbs/kraken2_standard_20221209",
        dir = "04-analysis/kraken2"
    threads: 16
    shell:
        """
         kraken2 --db {params.db} \
             --threads {threads} \
             --classified-out {params.dir}/{wildcards.sample}_classifiedreads.txt \
             --output - \
             --report {output} \
             --confidence 0.15 \
             --use-names \
             {input.fq} && gzip {params.dir}/{wildcards.sample}_classifiedreads.txt
        """
