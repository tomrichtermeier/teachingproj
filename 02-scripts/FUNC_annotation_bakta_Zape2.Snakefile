################################################################################
# Project: Training project
# Part: Functional profiling
# Step: Annotation of Zape2 contigs using Bakta
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#
# Alex Huebner, 12/04/23
################################################################################

from glob import glob
import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### ASSEMBLIES ################################################################
# MEGAHIT
MEGAHIT = {os.path.basename(fn).replace(".contigs.fa.gz", ""): fn
           for fn in glob("04-analysis/Zape2/assembly/Assembly/MEGAHIT/*.contigs.fa.gz")}
# SPADES
SPADES = {os.path.basename(fn).replace("_scaffolds.fasta.gz", ""): fn
          for fn in glob("04-analysis/Zape2/assembly/Assembly/SPAdes/*_scaffolds.fasta.gz")}
# Merge
SAMPLES = {**MEGAHIT, **SPADES}
################################################################################

rule all:
    input:
        expand("04-analysis/Zape2/bakta/{sample}.gff3", sample=SAMPLES)

rule download_bakta_db:
    # ln -s /mnt/archgen/microbiome_paleobiotech/calcBGCecoevo/03-data/refdbs/bakta 03-data/refdbs/ && touch 03-data/refdbs/bakta/downloaded
    # links and then creates dummy
    output:
        touch("03-data/refdbs/bakta/downloaded")
    message: "Download and prepare the reference DB for Bakta"
    conda: "ENVS_bakta.yaml"
    resources:
        mem = 4,
        cores = 1
    params:
        dir = "03-data/refdbs/bakta"
    shell:
        """
        bakta_db download --output {params.dir}
        """

rule bakta:
    # this rule uses singularity containers
    # give access to database because container doesnt have access to it in default   
    # executed with --use-singularity --singularity-args "-B /mnt/archgen/microbiome_paleobiotech/calcBGCecoevo/03-data/refdbs/bakta/db"
    input:
        lambda wildcards: SAMPLES[wildcards.sample],
        "03-data/refdbs/bakta/downloaded"
    output:
        "04-analysis/Zape2/bakta/{sample}.gff3"
    message: "Annotate contigs using BAKTA: {wildcards.sample}"
    container: "/mnt/archgen/tools/singularity/containers/depot.galaxyproject.org-singularity-bakta-1.7.0--pyhdfd78af_0.img"
    resources:
        mem = 128,
        cores = 8
    params:
        prefix = lambda wildcards: wildcards.sample,
        outdir = "04-analysis/Zape2/bakta",
        dbdir = "03-data/refdbs/bakta/db",
        extra = "--keep-contig-headers --meta --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-plot --skip-ori --skip-gap --skip-crispr"
    threads: 8
    wrapper:
        "file:///home/alexander_huebner/github/snakemake-wrappers/bio/bakta/bakta"

################################################################################

#### Compress files ############################################################

rule compress_files:
    params:
        dir = "04-analysis/Zape2/bakta"
    shell:
        "pigz -p {params.dir}/*"

################################################################################
