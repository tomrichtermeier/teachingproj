################################################################################
# Project: Training project
# Part: Functional profiling
# Step: Identify BGCs using antiSMASH for Zape2
#
# Due to difficulties of running antiSMASH from within nf-core/funcscan, we
# will run prodigal and antiSMASH outside of it.
#
# Dependent on:
#   - ASMB_nfcore_mag_Zape2.Snakefile
#
# Alex Huebner, 24/05/23
################################################################################

import os

import pyfastx
import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
MEGAHIT = {s: f"04-analysis/Zape2/assembly/Assembly/MEGAHIT/{s}.contigs.fa.gz"
           for s in glob_wildcards("04-analysis/Zape2/assembly/Assembly/MEGAHIT/{sample}.contigs.fa.gz")[0]}
SPADES = {s: f"04-analysis/Zape2/assembly/Assembly/SPAdes/{s}_contigs.fasta.gz"
          for s in glob_wildcards("04-analysis/Zape2/assembly/Assembly/SPAdes/{sample}_contigs.fasta.gz")[0]}
SAMPLES = {**MEGAHIT, **SPADES}
print(SAMPLES)
################################################################################

rule all:
    input:
        expand("04-analysis/Zape2/antismash/{sample}/{sample}.tsv", sample=SAMPLES)

rule decompress_fasta:
    output:
        pipe("tmp/pyrodigal/{sample}.fasta")
    message: "Decompress FastA: {wildcards.sample}"
    resources:
        mem = 2,
        cores = 1
    params:
        fasta = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        "gunzip -c {params.fasta} > {output}"

rule pyrodigal:
    input:
        "tmp/pyrodigal/{sample}.fasta"
    output:
        fna = "04-analysis/Zape2/pyrodigal/{sample}.fna",
        faa = "04-analysis/Zape2/pyrodigal/{sample}.faa",
        gff = "04-analysis/Zape2/pyrodigal/{sample}.gff",
        score = "04-analysis/Zape2/pyrodigal/{sample}.score.txt"
    message: "Identify CDS using pyrodigal: {wildcards.sample}"
    container: "/mnt/archgen/tools/singularity/containers/depot.galaxyproject.org-singularity-pyrodigal-2.1.0--py310h1425a21_0.img"
    resources:
        mem = 16,
        cores = 1
    shell:
        """
        pyrodigal -i {input} \
            -p meta \
            -o {output.gff} \
            -d {output.fna} \
            -a {output.faa} \
            -s {output.score}
        """
    
rule antismash:
    input:
        "04-analysis/Zape2/pyrodigal/{sample}.gff"
    output:
        zip = "04-analysis/Zape2/antismash/{sample}/{sample}.zip",
        gbk = "04-analysis/Zape2/antismash/{sample}/{sample}.gbk"
    message: "Run antiSMASH: {wildcards.sample}"
    #container: "docker://antismash/standalone"
    singularity: "/mnt/archgen/users/huebner/containers/antismash_standalone_v7.0.0.sif"
    resources:
        mem = 50,
        cores = 8
    params:
        fasta = lambda wildcards: SAMPLES[wildcards.sample],
        outdir = "04-analysis/Zape2/antismash/{sample}"
    threads: 8
    log: "04-analysis/Zape2/antismash/{sample}/antiSMASH.log"
    shell:
        """
        antismash \
            --clusterhmmer --tigrfam --asf --cc-mibig --cb-general --rre \
            -c {threads} \
            --genefinding-tool none \
            --genefinding-gff3 {input} \
            --cb-knownclusters \
            --output-dir {params.outdir} \
            --output-basename {wildcards.sample} \
            --minlength 3000 --allow-long-headers \
            --logfile {log} \
            {params.fasta} # && \
        """
        #find {params.outdir} -mindepth 1 -maxdepth 1 -type d -exec rm -r {{}} + && \
        #find . -type f ! \( -name "antiSMASH.log" -o -name "{wildcards.sample}.gbk" -o -name "{wildcards.sample}.zip" \) -exec rm {{}} +

rule comBGC:
    input:
        "04-analysis/Zape2/antismash/{sample}/{sample}.gbk"
    output:
        "04-analysis/Zape2/antismash/{sample}/{sample}.tsv"
    message: "Parse antiSMASH output to TSV: {wildcards.sample}"
    conda: "ENVS_comBGC.yaml"
    resources:
        mem = 4,
        cores = 1
    params:
        outdir = "04-analysis/Zape2/antismash/{sample}"
    threads: 1
    shell:
        """
        $HOME/.nextflow/assets/nf-core/funcscan/bin/comBGC.py -i {input} -o {params.outdir} && \
        mv {params.outdir}/combgc_summary.tsv {output}
        """
