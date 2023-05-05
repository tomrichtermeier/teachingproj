################################################################################
# Project: Training project
# Part: Phylogenetic analyses
# Step: Maximum likelihood tree of P. copri
#
# Dependent on:
#   - REFG_dRep_Pcopri.Snakefile
#
# Alex Huebner, 03/05/23
################################################################################

import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### GENOMES ###################################################################
# Reconstructed genomes
RECON_GENOMES, = glob_wildcards("04-analysis/refgenome_reconst/{genome}.freebayes_conserv.fasta")
RECON_GENOMES = [g for g in RECON_GENOMES if g[:3] == "SRS" or g[:6] == "ERS418"]
# Comparative genomes
SGB_GENOMES, = glob_wildcards("03-data/refgenomes/pcopri/sgb/{sgb}.fa")
REFSEQ_GENOMES, = glob_wildcards("03-data/refgenomes/pcopri/refseq/{refseq}.fa")
GENOMES = RECON_GENOMES + SGB_GENOMES + REFSEQ_GENOMES
################################################################################

#### Auxilliary functions ######################################################

def return_fa_path(wildcards):
    if wildcards.genome[:3] == "SRS" or wildcards.genome[:6] == "ERS418":
        return f"{os.getcwd()}/04-analysis/refgenome_reconst/{wildcards.genome}.freebayes_conserv.fasta"
    elif wildcards.genome[:3] == "GCA":
        return f"{os.getcwd()}/03-data/refgenomes/pcopri/refseq/{wildcards.genome}.fa"
    else:
        return f"{os.getcwd()}/03-data/refgenomes/pcopri/sgb/{wildcards.genome}.fa"

################################################################################

localrules: link_fasta

rule all:
    input:
        "04-analysis/phylophlan3/pcopri/RAxML_bestTree.genomes_refined.tre"

#### Prepare genomes ###########################################################

rule link_fasta:
    output:
        "tmp/phylophlan_pcopri/genomes/{genome}.fa"
    message: "Link FastA file into temp folder: {wildcards.genome}"
    params:
        fa = lambda wildcards: return_fa_path(wildcards)
    threads: 1
    shell:
        "ln -s {params.fa} {output}"

################################################################################

#### Prepare marker-gene database ##############################################

rule setup_database:
    output:
        "tmp/phylophlan_pcopri/marker_genes/s__Prevotella_copri/s__Prevotella_copri.fna"
    message: "Construct marker-gene databae for P. copri from UniProt"
    conda: "ENVS_phylophlan.yaml"
    resources:
        mem = 4,
        cores =1
    params:
        dir = "tmp/phylophlan_pcopri/marker_genes"
    shell:
        """
        phylophlan_setup_database \
            --get_core_proteins s__Prevotella_copri \
            -o {params.dir} \
            -t n \
            --verbose
        """

################################################################################

#### Run PhyloPhlan ############################################################

rule write_config_nucleotide:
    output:
        "04-analysis/phylophlan3/pcopri/nucleotide_config.txt"
    message: "Write config file for nucleotide analysis for PhyloPhlAn3"
    conda: "ENVS_phylophlan.yaml"
    shell:
        """
        phylophlan_write_config_file \
            -o {output} \
            -d n\
            --db_dna makeblastdb \
            --map_dna blastn \
            --msa mafft \
            --trim trimal \
            --tree1 fasttree \
            --tree2 raxml
    """

rule phylophlan:
    input:
        configfn = "04-analysis/phylophlan3/pcopri/nucleotide_config.txt",
        db = "tmp/phylophlan_pcopri/marker_genes/s__Prevotella_copri/s__Prevotella_copri.fna",
        samples = expand("tmp/phylophlan_pcopri/genomes/{genome}.fa", genome=GENOMES)
    output:
        "04-analysis/phylophlan3/pcopri/RAxML_bestTree.genomes_refined.tre"
    message: "Run PhyloPhlAn3 analysis for species P. copri"
    conda: "ENVS_phylophlan.yaml"
    resources:
        mem = 72,
        cores = 16
    params:
        inputdir = "tmp/phylophlan_pcopri/genomes",
        outputdir = "04-analysis/phylophlan3/pcopri",
        dbdir = "tmp/phylophlan_pcopri/marker_genes"
    log: "04-analysis/phylophlan3/pcopri/phylophlan.log"
    threads: 16
    shell:
        """
        phylophlan \
            -i {params.inputdir} \
            -o {params.outputdir} \
            -d s__Prevotella_copri \
            --databases_folder {params.dbdir} \
            -t n \
            -f {input.configfn} \
            --genome_extension .fa \
            --nproc {threads} \
            --diversity low \
            --fast \
            --verbose 2>&1 | tee {log}
        """

################################################################################
