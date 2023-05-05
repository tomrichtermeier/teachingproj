################################################################################
# Project: Training project
# Part: Preparation of the data
# Step: Run dRep on T. succinifaciens genomes
################################################################################

from glob import glob
import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### COMPARATIVE DATA ##########################################################
SGB_OVERVIEW = pd.read_csv("01-documentation/41587_2023_1688_MOESM3_ESM_TableS1.tsv", sep="\t",
                           usecols=['Label', 'ID', 'SGB centroid', 'Assigned taxonomy']) \
    .query('Label == "SGB"') \
    .query('ID.isin([3539, 3542, 3546])')
SGBS = {'RampelliS_2015__H10__bin.14': 3539,
        'CM_madagascar__V07_01_2FE__bin.12': 3539,
        'BritoIL_2016__M2.38.ST__bin.21': 3539,
        'ZellerG_2014__CCIS35100175ST-4-0__bin.37': 3542,
        'SmitsSA_2017__TZ_65642__bin.4': 3542,
        'FengQ_2015__SID31160__bin.51': 3542,
        'SchirmerM_2016__G89199__bin.16': 3546,
        'Obregon-TitoAJ_2015__SM40__bin.39': 3546,
        'LeChatelierE_2013__MH0189__bin.74': 3546,
        'BritoIL_2016__M2.40.ST__bin.49': 3546}
NCBI_REFSEQ = ['GCA_000195275']
################################################################################

#### Auxilliary functions ######################################################

NCBI_REFSEQ_URLS = {'GCA_000195275': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/275/GCF_000195275.1_ASM19527v1/GCF_000195275.1_ASM19527v1_genomic.fna.gz'}

def return_url(wildcards):
        sgbid = SGBS[wildcards.genome]
        return f"https://s3-eu-west-1.amazonaws.com/sgb-genomes/data/genomes/{sgbid}/{wildcards.genome}.fa"

################################################################################


rule all:
    input:
        "05-results/REFG_dRep_Tsuccinifaciens_clusters.tsv"

#### Download comparative data #################################################

rule download_from_aws:
    output:
        "03-data/refgenomes/tsuccinifaciens/sgb/{genome}.fa"
    message: "Download SGB centroid: {wildcards.genome}"
    params:
        url = lambda wildcards: return_url(wildcards)
    shell:
        "wget -O {output} {params.url}"

rule download_from_refseq:
    output:
        "03-data/refgenomes/tsuccinifaciens/refseq/{genome}.fa"
    message: "Download genome from RefSeq: {wildcards.genome}"
    params:
        url = lambda wildcards: NCBI_REFSEQ_URLS[wildcards.genome],
        fna = "03-data/refgenomes/tsuccinifaciens/refseq/{genome}.fna.gz"
    shell:
        "wget -O {params.fna} {params.url} && gunzip -c {params.fna} > {output}"


################################################################################

#### dRep ######################################################################

rule prepare_genome_list:
    input:
        sgb = expand("03-data/refgenomes/tsuccinifaciens/sgb/{genome}.fa", genome=SGBS),
        refseq = expand("03-data/refgenomes/tsuccinifaciens/refseq/{genome}.fa", genome=NCBI_REFSEQ)
    output:
        "04-analysis/dRep/tsuccinifaciens_genomelist.txt"
    message: "Write list of T. succinifaciens genomes to be analysed by dRep"
    params:
        dir = "04-analysis/refgenome_reconst"
    run:
        reconst_genome_fns = [f for f in glob(f"{params.dir}/*.fasta")
                              if os.path.basename(f)[:3] == "SRS" or os.path.basename(f)[:6] == "ERS418"]
        with open(output[0], "wt") as outfile:
            for g in reconst_genome_fns:
                outfile.write(f"{os.getcwd()}/{g}\n")
            for g in input.sgb:
                outfile.write(f"{os.getcwd()}/{g}\n")
            for g in input.refseq:
                outfile.write(f"{os.getcwd()}/{g}\n")

rule dRep:
    input:
        "04-analysis/dRep/tsuccinifaciens_genomelist.txt"
    output:
        "04-analysis/dRep/tsuccinifaciens/data_tables/Widb.csv"
    message: "Run dRep on T. succinifaciens genomes"
    conda: "ENVS_dRep.yaml"
    resources:
        mem = 72,
        cores = 4
    params:
        outprefix = "04-analysis/dRep/tsuccinifaciens"
    threads: 8
    shell:
        """
        dRep dereplicate -p {threads} {params.outprefix} -g {input} -comp 50 --S_algorithm ANImf -pa 0.95 -sa 0.97
        """

rule export_table:
    input:
        "04-analysis/dRep/tsuccinifaciens/data_tables/Widb.csv"
    output:
        "05-results/REFG_dRep_Tsuccinifaciens_clusters.tsv"
    message: "Export final cluster table"
    run:
        pd.read_csv(input[0], sep=",") \
            .to_csv(output[0], sep="\t", index=False, float_format="%.2f")

################################################################################
