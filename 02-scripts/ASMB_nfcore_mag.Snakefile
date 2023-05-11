################################################################################
# Project: Training project
# Part: De novo assembly and genome binning
# Step: De novo assembly for all downloaded data with nf-core/mag
#
# Dependent on:
#   - PREP_remove_hostDNA_Zape2.Snakefile
#
# Alex Huebner, 09/05/23
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

localrules: prepare_samplesheet

rule all:
    input:
        "04-analysis/nfcore_mag/nfcore_mag.done"

rule prepare_samplesheet:
    output:
        "04-analysis/nfcore_mag/nfcore_mag_samplesheet.csv"
    message: "Create the samplesheet for nf-core/mag"
    threads: 1
    run:
        tbl = pd.DataFrame.from_dict({'sample': SAMPLES})
        tbl['group'] = list(range(6))
        tbl['short_reads_1'] = [f"03-data/processed_data/{sample}_1.fastq.gz" for sample in SAMPLES]
        tbl['short_reads_2'] = [f"03-data/processed_data/{sample}_2.fastq.gz" for sample in SAMPLES]
        tbl['long_reads'] = ""
        tbl.to_csv(output[0], sep=",", index=False)

rule nfcore_mag:
    input:
        "04-analysis/nfcore_mag/nfcore_mag_samplesheet.csv"
    output:
        touch("04-analysis/nfcore_mag/nfcore_mag.done")
    message: "Assemble samples with nf-core/mag"
    params:
        outdir = "04-analysis/nfcore_mag/assembly"
    shell:
        """
        nextflow run nf-core/mag -r 2.3.0 \
            -profile eva,archgen \
            --input {input} \
            --outdir {params.outdir} \
            --skip_clipping \
            --skip_prodigal \
            --binning_map_mode own \
            --min_contig_size 1000 \
            --binqc_tool checkm \
            --save_checkm_data \
            --refine_bins_dastool \
            --postbinning_input both \
            --run_gunc \
            --gunc_save_db \
            --ancient_dna
        """
