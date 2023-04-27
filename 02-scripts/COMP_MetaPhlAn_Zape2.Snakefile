################################################################################
# Project: Training project
# Part: Composition analysis
# Step: Taxonomic profiling of Zape2 with MetaPhlAn
#
# Dependent on:
#   - PREP_remove_hostDNA.Snakefile
#
# Alex Huebner, 19/04/23
################################################################################

import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

rule all:
    input:
        "04-analysis/metaphlan4/Zape2.metaphlan.profile.txt"
    
rule bam2fq:
    output:
        pipe("tmp/metaphlan/Zape2.fastq")
    message: "Extract unmapped reads and convert into FastQ"
    conda: "ENVS_samtools.yaml"
    params:
        bam = "../week2/eager_Zape2/mapping/bwa/ERR3678612_PE.mapped.bam"
    shell:
        """
        samtools view -uh -f 4 {params.bam} | samtools fastq - > {output}
        """

rule metaphlan4_install_db:
    output:
        touch("04-analysis/metaphlan4/installed_database.done")
    message: "Install the CHOCOPhlAn database"
    conda: "ENVS_MetaPhlAn4.yaml"
    shell:
        "metaphlan --install"
    
rule metaphlan:
    input:
        db = "04-analysis/metaphlan4/installed_database.done",
        fq = "tmp/metaphlan/Zape2.fastq"
    output:
        sam = "04-analysis/metaphlan4/Zape2.metaphlan.sam.bz2",
        profile = "04-analysis/metaphlan4/Zape2.metaphlan.profile.txt"
    message: "Run MetaPhlAn4 with default settings for sample Zape2"
    conda: "ENVS_MetaPhlAn4.yaml"
    resources:
        mem = 24,
        cores = 8
    params:
        tmp_dir = "tmp/metaphlan"
    threads: 8
    shell:
        """
        metaphlan \
            {input.fq} \
            --input_type fastq \
            --tmp_dir {params.tmp_dir} \
            --force \
            --index mpa_vOct22_CHOCOPhlAnSGB_202212 \
            --ignore_eukaryotes \
            -t rel_ab_w_read_stats \
            --sample_id Zape2 \
            --read_min_len 35 \
            -s {output.sam} \
            -o {output.profile} \
            --nproc {threads}
        """
