################################################################################
# Project: Training project
# Part: Quality evaluation
# Step: Determine the amount of ancient DNA damage
#
# Dependent on:
#   - PREP_remove_hostDNA.Snakefile
#   - COMP_Kraken2_Bracken.Snakefile
#
# Alex Huebner, 25/04/23
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

GENOME_URLS = {'Pcopri': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/735/445/GCF_020735445.1_ASM2073544v1/GCF_020735445.1_ASM2073544v1_genomic.fna.gz',
               'Bsubtilis': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_cds_from_genomic.fna.gz',
               'Smaltophilia': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/475/405/GCF_900475405.1_44087_C01/GCF_900475405.1_44087_C01_genomic.fna.gz'}

wildcard_constraints:
    sample = "[ES]RS[0-9]+"

rule all:
    input:
        "05-results/QUAL_damageprofile_Pcopri_Smaltophilia.tsv"

#### Prepare sequencing data ###################################################

rule decompress_fasta:
    output:
        temp("tmp/damageprofiles/{genome}.fna")
    message: "Decompress the FastA file: {wildcards.genome}"
    params:
        url = lambda wildcards: GENOME_URLS[wildcards.genome]
    shell:
        "wget -O - {params.url} | gunzip > {output}"

rule bgzip_tabix:
    input:
        "tmp/damageprofiles/{genome}.fna"
    output:
        fasta = "tmp/damageprofiles/{genome}.fna.gz",
        fai = "tmp/damageprofiles/{genome}.fna.gz.fai"
    message: "Compress the FastA file: {wildcards.genome}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 4,
        cores = 1
    threads: 1
    shell:
        """
        bgzip {input}
        samtools faidx {output.fasta}
        """

rule bowtie2_index:
    input:
        fasta = "tmp/damageprofiles/{genome}.fna.gz",
        fai = "tmp/damageprofiles/{genome}.fna.gz.fai"
    output:
        "tmp/damageprofiles/{genome}.1.bt2"
    message: "Index for BowTie2 alignment: {wildcards.genome}"
    conda: "ENVS_bowtie2.yaml"
    resources:
        mem = 4,
        cores = 1
    params:
        prefix = "tmp/damageprofiles/{genome}"
    threads: 1
    shell:
        """
        bowtie2-build --threads {threads} \
                {input.fasta} {params.prefix}
        """

################################################################################

#### Align data ################################################################

rule bowtie2:
    input:
        lambda wildcards: "tmp/damageprofiles/Pcopri.1.bt2" if wildcards.sample[:3] == "SRS" or wildcards.sample[:6] == "ERS418" else "tmp/damageprofiles/Smaltophilia.1.bt2"
    output:
        pipe("tmp/damageprofiles/{sample}.sam")
    message: "Align sequences against reference genomes using BowTie2: {wildcards.sample}"
    conda: "ENVS_bowtie2.yaml"
    resources:
        mem = 8,
        cores = 8
    params:
        index = lambda wildcards: "tmp/damageprofiles/Pcopri" if wildcards.sample[:3] == "SRS" or wildcards.sample[:6] == "ERS418" else "tmp/damageprofiles/Smaltophilia",
        pe1 = "03-data/processed_data/{sample}_1.fastq.gz",
        pe2 = "03-data/processed_data/{sample}_2.fastq.gz"
    threads: 8
    shell:
        """
        bowtie2 -p {threads} --very-sensitive -x {params.index} \
            -1 {params.pe1} -2 {params.pe2} > {output}
        """

rule sam2bam:
    input:
        "tmp/damageprofiles/{sample}.sam"
    output:
        pipe("tmp/damageprofiles/{sample}.bam")
    message: "Convert SAM to BAM format: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 2,
        cores = 0
    shell:
        "samtools view -Su -e '!flag.unmap || !flag.munmap' {input} > {output}"

rule samtools_fixmate:
    input:
        "tmp/damageprofiles/{sample}.bam"
    output:
        pipe("tmp/damageprofiles/{sample}.fixmate.bam")
    message: "Fix mate flags: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 2,
        cores = 0
    shell:
        "samtools fixmate -mu {input} {output}"

rule samtools_sort:
    input:
        "tmp/damageprofiles/{sample}.fixmate.bam"
    output:
        pipe("tmp/damageprofiles/{sample}.sorted.bam")
    message: "Sort BAM file by coordinate: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 2,
        cores = 0
    shell:
        "samtools sort -u -o {output} {input}"

rule samtools_calmd:
    input:
        "tmp/damageprofiles/{sample}.sorted.bam"
    output:
        "04-analysis/damageprofiles/{sample}.calmd.bam"
    message: "Calculate the MD tag: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 1
    params:
        fa = lambda wildcards: "tmp/damageprofiles/Pcopri.fna.gz" if wildcards.sample[:3] == "SRS" or wildcards.sample[:6] == "ERS418" else "tmp/damageprofiles/Smaltophilia.fna.gz"
    shell:
        "samtools calmd -b {input} {params.fa} > {output}"

rule samtools_index:
    input:
        "04-analysis/damageprofiles/{sample}.calmd.bam"
    output:
        "04-analysis/damageprofiles/{sample}.calmd.bam.bai"
    message: "Index the BAM file: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 4,
        cores = 1
    shell:
        "samtools index {input}"

################################################################################

#### Quantify damage ###########################################################

rule damageprofiler:
    input:
        bam = "04-analysis/damageprofiles/{sample}.calmd.bam",
        bai = "04-analysis/damageprofiles/{sample}.calmd.bam.bai"
    output:
        "04-analysis/damageprofiles/{sample}/5p_freq_misincorporations.txt"
    message: "Profile the aDNA damage using damageprofiler: {wildcards.sample}"
    conda: "ENVS_damageprofiler.yaml"
    resources:
        mem = 8
    params:
        outdir = "04-analysis/damageprofiles/{sample}"
    shell:
        """
        damageprofiler -i {input.bam} \
            -o {params.outdir}
        """

rule summarise_damageprofiler:
    input:
        expand("04-analysis/damageprofiles/{sample}/5p_freq_misincorporations.txt", sample=SAMPLES)
    output:
        "05-results/QUAL_damageprofile_Pcopri_Smaltophilia.tsv"
    message: "Summarise the substitution frequency at the 5' end"
    run:
        damage = pd.concat([pd.read_csv(fn, sep="\t", skiprows=3) \
                                .assign(sample=os.path.basename(os.path.dirname(fn)))
                            for fn in input])
        damage['genome'] = ["Pcopri" if s[:3] == "SRS" or s[:6] == "ERS418" else "Smaltophilia" for s in damage['sample'].values]
        damage.iloc[:, [-2, -1] + list(range(0, 13))] \
            .to_csv(output[0], sep="\t", index=False, float_format="%.4f")

################################################################################
