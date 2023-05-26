################################################################################
# Project: Training project
# Part: Functional profiling
# Step: antiSMASH
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#
# Alex Huebner, 26/05/23
################################################################################

import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
METASPADES = {s: f"{os.getcwd()}/04-analysis/metagenome_assembly/alignment/metaspades/{s}-metaspades.fasta.gz"
              for s in glob_wildcards("04-analysis/metagenome_assembly/alignment/metaspades/{sample}-metaspades.fasta.gz")[0]}
SAMPLES = METASPADES
################################################################################

rule all:
    input:
        expand("04-analysis/funcscan/results/bgc/antismash/{sample}/{sample}.tsv", sample=SAMPLES)

rule decompress_fasta:
    output:
        pipe("04-analysis/funcscan/results/annotation/pyrodigal/{sample}.fasta")
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
        "04-analysis/funcscan/results/annotation/pyrodigal/{sample}.fasta"
    output:
        fna = "04-analysis/funcscan/results/annotation/pyrodigal/{sample}.fna",
        faa = "04-analysis/funcscan/results/annotation/pyrodigal/{sample}.faa",
        gff = "04-analysis/funcscan/results/annotation/pyrodigal/{sample}.gff"
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
            -a {output.faa}
        """
    
rule antismash:
    input:
        "04-analysis/funcscan/results/annotation/pyrodigal/{sample}.gff"
    output:
        zip = "04-analysis/funcscan/results/bgc/antismash/{sample}/{sample}.zip",
        gbk = "04-analysis/funcscan/results/bgc/antismash/{sample}/{sample}.gbk"
    message: "Run antiSMASH: {wildcards.sample}"
    #container: "docker://antismash/standalone"
    singularity: "/mnt/archgen/users/huebner/containers/antismash_standalone_v7.0.0.sif"
    resources:
        mem = 50,
        cores = 8
    params:
        fasta = lambda wildcards: SAMPLES[wildcards.sample],
        outdir = "04-analysis/funcscan/results/bgc/antismash/{sample}"
    threads: 8
    log: "04-analysis/funcscan/results/bgc/antismash/{sample}/antiSMASH.log"
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
        "04-analysis/funcscan/results/bgc/antismash/{sample}/{sample}.gbk"
    output:
        "04-analysis/funcscan/results/bgc/antismash/{sample}/{sample}.tsv"
    message: "Parse antiSMASH output to TSV: {wildcards.sample}"
    conda: "ENVS_comBGC.yaml"
    resources:
        mem = 4,
        cores = 1
    params:
        outdir = "04-analysis/funcscan/results/bgc/antismash/{sample}"
    threads: 1
    shell:
        """
        $HOME/.nextflow/assets/nf-core/funcscan/bin/comBGC.py -i {input} -o {params.outdir} && \
        mv {params.outdir}/combgc_summary.tsv {output}
        """
