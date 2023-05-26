################################################################################
# Project: Training project
# Part: Functional profiling
# Step: Annotation using FuncScan
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#
# Alex Huebner, 25/05/23
################################################################################

import os

import pandas as pd

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
        "04-analysis/funcscan/funcscan_pyrodigal.done"

rule create_samplesheet:
    output:
        "04-analysis/funcscan/samplesheet.csv"
    message: "Create the samplesheet for Funcscan"
    resources:
        mem = 4,
        cores = 1
    threads: 1
    run:
        pd.DataFrame.from_dict(SAMPLES, orient="index", columns=['fasta']) \
            .to_csv(output[0], sep=",", index=True, index_label="sample")

rule download_deepbgc_db:
    output:
        touch("03-data/refdbs/deepbgc/downloaded")
    message: "Download and prepare the reference DB for DeepBGC"
    conda: "ENVS_deepbgc.yaml"
    resources:
        mem = 4,
        cores = 1
    params:
        dir = "03-data/refdbs/deepbgc"
    shell:
        """
        DEEPBGC_DOWNLOADS_DIR={params.dir} deepbgc download
        """

rule funcscan_pyrodigal:
    input:
        samples = "04-analysis/funcscan/samplesheet.csv",
        deepbgc = "03-data/refdbs/deepbgc/downloaded"
    output:
        touch("04-analysis/funcscan/funcscan_pyrodigal.done")
    message: "Annotate contigs using nf-core/funcscan with pyrodigal"
    params:
        outdir = "04-analysis/funcscan/results",
        deepbgc_db = f"{os.getcwd()}/03-data/refdbs/deepbgc/",
        custom_config = "01-documentation/funcscan_contigs.conf"
    shell:
        """
        nextflow run nf-core/funcscan -r 1.1.1 \
            -profile eva,archgen \
            --input {input.samples} \
            --outdir {params.outdir} \
            --run_bgc_screening \
            --annotation_tool pyrodigal --save_annotations \
            --bgc_skip_antismash \
            --bgc_deepbgc_database {params.deepbgc_db} \
            --bgc_skip_hmmsearch \
            -c {params.custom_config}
        """
