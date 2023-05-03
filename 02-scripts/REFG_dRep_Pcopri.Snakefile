################################################################################
# Project: Training project
# Part: Preparation of the data
# Step: Run dRep on Prevotella copri genomes
#
# Dependent on:
#   - REFG_genome_reconstruction_Pcopri_Smaltophilia.Snakefile
#
# Alex Huebner, 02/05/23
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
    .query('ID > 1610 and ID <= 1656') #Filtering IDs that are interesting for Pcopri
# Replace non-existant SGBs, Problem: not all centroids are publuished on same server, has to be picked out manually (picked from other server)
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1612, 'SGB centroid'] = "LiJ_2014__V1.CD6-0-PT__bin.41"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1615, 'SGB centroid'] = "VincentC_2016__MM019.3__bin.15"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1617, 'SGB centroid'] = "NielsenHB_2014__O2_UC49_0__bin.6"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1623, 'SGB centroid'] = "ChengpingW_2017__AS133raw__bin.39"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1624, 'SGB centroid'] = "CM_madagascar__A29_01_1FE__bin.15"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1626, 'SGB centroid'] = "ZeeviD_2015__PNP_Validation_85__bin.4"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1636, 'SGB centroid'] = "CM_madagascar__V53_02_1FE__bin.39"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1639, 'SGB centroid'] = "QinN_2014__LD-75__bin.24"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1642, 'SGB centroid'] = "FengQ_2015__SID530743__bin.21"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1644, 'SGB centroid'] = "ZeeviD_2015__PNP_Main_200__bin.56"
SGB_OVERVIEW.loc[SGB_OVERVIEW['ID'] == 1653, 'SGB centroid'] = "LiJ_2017__H1M513807__bin.23"
SGB_OVERVIEW = SGB_OVERVIEW.set_index('SGB centroid')
# List genomes
SGBS = SGB_OVERVIEW.index.tolist()
################################################################################



#### Auxilliary functions ######################################################

NCBI_REFSEQ_URLS = {'GCA_000431915': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/431/915/GCF_000431915.1_MGS386/GCF_000431915.1_MGS386_genomic.fna.gz',
                    'GCA_900555035': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/555/035/GCA_900555035.1_UMGS1787/GCA_900555035.1_UMGS1787_genomic.fna.gz',
                    'GCA_900551985': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/551/985/GCA_900551985.1_UMGS1464/GCA_900551985.1_UMGS1464_genomic.fna.gz'}

def return_url(wildcards):
        sgbid = SGB_OVERVIEW.at[wildcards.genome, 'ID']
        return f"https://s3-eu-west-1.amazonaws.com/sgb-genomes/data/genomes/{sgbid}/{wildcards.genome}.fa" #Data are stored on cloud server 

################################################################################


rule all:
    input:
        "05-results/REFG_dRep_Pcopri_clusters.tsv"

#### Download comparative data #################################################

rule download_from_aws:
    output:
        "03-data/refgenomes/pcopri/sgb/{genome}.fa"
    message: "Download SGB centroid: {wildcards.genome}"
    params:
        url = lambda wildcards: return_url(wildcards)
    shell:
        "wget -O {output} {params.url}"

rule download_from_refseq:
    output:
        "03-data/refgenomes/pcopri/refseq/{genome}.fa"
    message: "Download genome from RefSeq: {wildcards.genome}"
    params:
        url = lambda wildcards: NCBI_REFSEQ_URLS[wildcards.genome],
        fna = "03-data/refgenomes/pcopri/refseq/{genome}.fna.gz"
    shell:
        "wget -O {params.fna} {params.url} && gunzip -c {params.fna} > {output}"


################################################################################

#### dRep ######################################################################

rule prepare_genome_list:
    input:
        sgb = expand("03-data/refgenomes/pcopri/sgb/{genome}.fa", genome=[s for s in SGBS if not s.startswith("GCA")]),
        refseq = expand("03-data/refgenomes/pcopri/refseq/{genome}.fa", genome=[s for s in SGBS if s.startswith("GCA")])
    output:
        "04-analysis/dRep/pcopri_genomelist.txt"
    message: "Write list of P. copri genomes to be analysed by dRep"
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
        "04-analysis/dRep/pcopri_genomelist.txt"
    output:
        "04-analysis/dRep/pcopri/data_tables/Widb.csv"
    message: "Run dRep on P. copri genomes"
    conda: "ENVS_dRep.yaml"
    resources:
        mem = 72,
        cores = 4
    params:
        outprefix = "04-analysis/dRep/pcopri"
    threads: 8
    shell:
        """
        dRep dereplicate -p {threads} {params.outprefix} -g {input} -comp 50 --S_algorithm ANImf -pa 0.95 -sa 0.99
        """

rule export_table:
    input:
        "04-analysis/dRep/pcopri/data_tables/Widb.csv"
    output:
        "05-results/REFG_dRep_Pcopri_clusters.tsv"
    message: "Export final cluster table"
    run:
        pd.read_csv(input[0], sep=",") \
            .to_csv(output[0], sep="\t", index=False, float_format="%.2f")

################################################################################
