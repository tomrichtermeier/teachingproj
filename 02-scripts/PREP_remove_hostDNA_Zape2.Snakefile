################################################################################
# Project: Training project
# Part: Preparation of the data
# Step: Extract non-host DNA of Zape2 using nf-core/eager
################################################################################

import os
import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES = {'Zape2_MinE2X': ['ERR3678612', 'ERR3678626'],
           'Zape2_PCMinE': ['ERR3678620', 'ERR3678629']}
################################################################################

#### Auxilliary functions ######################################################

def return_samples(wildcards):
    eager_tbl = pd.read_csv(checkpoints.amdirt_convert_eager.get(**wildcards).output[0], sep="\t")
    if wildcards.sample in SAMPLES:
        run_accs = SAMPLES[wildcards.sample]
    else:
        run_accs = eager_tbl['Library_ID'].tolist()
    return [f"tmp/eager_extract_unmapped/{racc}_{i}.fastq.gz" for racc in run_accs for i in range(3)]


def return_fastq_path(wildcards, mate):
    eager_tbl = pd.read_csv(checkpoints.amdirt_convert_eager.get(**wildcards).output[0], sep="\t")
    if wildcards.sample in SAMPLES:
        run_accs = SAMPLES[wildcards.sample]
    else:
        run_accs = eager_tbl['Library_ID'].tolist()
    return [f"tmp/eager_extract_unmapped/{racc}_{mate}.fastq.gz" for racc in run_accs]

################################################################################

rule all:
    input:
        "05-results/PREP_Nextflow_EAGER_Zape2_noReads.tsv"

rule subset_sample_table:
    output:
        "04-analysis/Zape2/amdirt_samplesheet.tsv"
    message: "Filter AncientMetagenomeDir for sample Zape2 from Hagan2019"
    params:
        url = 'https://raw.githubusercontent.com/SPAAM-community/AncientMetagenomeDir/master/ancientmetagenome-hostassociated/samples/ancientmetagenome-hostassociated_samples.tsv'
    run:
        pd.read_csv(params.url, sep="\t") \
            .query('sample_name == "Zape2" and project_name == "Hagan2019"') \
            .to_csv(output[0], sep="\t", index=False)

rule amdirt_convert_aspera:
    input:
        "04-analysis/Zape2/amdirt_samplesheet.tsv"
    output:
        "04-analysis/Zape2/AncientMetagenomeDir_aspera_download_script.sh"
    message: "Prepare a samplesheet for running Zape2 using nf-core/eager"
    params:
        outdir = "04-analysis/Zape2"
    shell:
        """
        AMDirT convert {input} ancientmetagenome-hostassociated \
            --aspera --output {params.outdir}
        """

rule aspera_download:
    input:
        "04-analysis/Zape2/AncientMetagenomeDir_aspera_download_script.sh"
    output:
        touch("04-analysis/Zape2/aspera_download.done")
    message: "Download Zape2 sequencing data with aspera"
    params:
        dir = "03-data/raw_data/Zape2",
        aspera_path = "/usr/local64/opt/aspera/connect"
    shell:
        """
        mkdir -p {params.dir}
        cd {params.dir}
        export ASPERA_PATH="{params.aspera_path}"
        bash ../../../{input}
        """

checkpoint amdirt_convert_eager:
    input:
        tsv = "04-analysis/Zape2/amdirt_samplesheet.tsv",
        aspera = "04-analysis/Zape2/aspera_download.done"
    output:
        "04-analysis/Zape2/AncientMetagenomeDir_nf_core_eager_input_table.tsv"
    message: "Prepare a samplesheet for running Zape2 using nf-core/eager"
    params:
        outdir = "04-analysis/Zape2",
        seqdatadir = "03-data\/raw_data\/Zape2"
    shell:
        """
        AMDirT convert {input.tsv} ancientmetagenome-hostassociated \
            --eager --output {params.outdir} && \
        sed -i "s/ERR/{params.seqdatadir}\/ERR/g" {output} && \
        sed -i "s/{params.seqdatadir}\/ERR/ERR/" {output}
        """

rule nfcore_eager:
    input:
        "04-analysis/Zape2/AncientMetagenomeDir_nf_core_eager_input_table.tsv"
    output:
        touch("04-analysis/Zape2/nfcore_eager.done")
    message: "Run nf-core/eager on dog palaeofaeces sample: Zape2"
    params:
    shell:
        """
        nextflow run nf-core/eager -r 2.4.6 \
          -profile eva,archgen \
          --input {input} \
          --fasta '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa' \
          --fasta_index '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa.fai' \
          --bwa_index '/mnt/archgen/Reference_Genomes/Human/hs37d5' \
          --seq_dict '/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.dict' \
          --skip_preseq \
          --skip_deduplication \
          --skip_damage_calculation \
          --skip_qualimap \
          --complexity_filter_poly_g \
          --skip_collapse \
          --bwaalno 1 --bwaalnl 32 \
          --outdir "04-analysis/Zape2/eager"
        """

rule samtools_sort_by_name:
    input:
        "04-analysis/Zape2/nfcore_eager.done"
    output:
        pipe("tmp/eager_extract_unmapped/{lib}.nsorted.bam")
    message: "Sort the BAM file by name: {wildcards.lib}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 12,
        cores = 4
    params:
        bam = "04-analysis/Zape2/eager/mapping/bwa/{lib}_PE.mapped.bam"
    threads: 4
    shell:
        """
        samtools sort -n -o {output} {params.bam}
        """
        
rule samtools_fixmate:
    input:
        "tmp/eager_extract_unmapped/{lib}.nsorted.bam"
    output:
        pipe("tmp/eager_extract_unmapped/{lib}.fixmate.bam")
    message: "Fix mate flags: {wildcards.lib}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 4
    threads: 4
    shell:
        """
        samtools fixmate -mcu -@ {threads} {input} {output}
        """

rule extract_unmapped_reads:
    input:
        "tmp/eager_extract_unmapped/{lib}.fixmate.bam"
    output:
        pe1 = temp("tmp/eager_extract_unmapped/{lib}_1.fastq.gz"),
        pe2 = temp("tmp/eager_extract_unmapped/{lib}_2.fastq.gz"),
        pe0 = temp("tmp/eager_extract_unmapped/{lib}_0.fastq.gz")
    message: "Extract all reads for which are not aligned in a proper pair and convert to fastq: {wildcards.lib}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 2
    threads: 2
    shell:
        """
        samtools view -uh -e '(flag.paired && (flag.unmap || flag.munmap)) || (!flag.paired && flag.unmap)' {input} | \
        samtools fastq -1 {output.pe1} \
                       -2 {output.pe2} \
                       -0 {output.pe0} -
        """
        
rule concat_fastqs:
    input:
        fqs = lambda wildcards: return_samples(wildcards)
    output:
        pe1 = "03-data/processed_data/Zape2/{sample}_1.fastq.gz",
        pe2 = "03-data/processed_data/Zape2/{sample}_2.fastq.gz",
        pe0 = "03-data/processed_data/Zape2/{sample}_0.fastq.gz"
    message: "Concatenate the FastQ files"
    params:
        tmpdir = "tmp/eager_extract_unmapped",
        outdir = "03-data/processed_data/Zape2",
        pe1 = lambda wildcards: " ".join(return_fastq_path(wildcards, 1)),
        pe2 = lambda wildcards: " ".join(return_fastq_path(wildcards, 2)),
        pe0 = lambda wildcards: " ".join(return_fastq_path(wildcards, 0))
    shell:
        """
        cat {params.pe1} > {output.pe1}
        cat {params.pe2} > {output.pe2}
        cat {params.pe0} > {output.pe0}
        """

rule count_reads:
    input:
        pe1 = "03-data/processed_data/Zape2/{sample}_1.fastq.gz",
        pe2 = "03-data/processed_data/Zape2/{sample}_2.fastq.gz",
        pe0 = "03-data/processed_data/Zape2/{sample}_0.fastq.gz"
    output:
        temp("03-data/processed_data/Zape2/{sample}.n")
    message: "Count the number of reads: {wildcards.sample}"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 2
    shell:
        """
        reads_PE1=$(bioawk -c fastx 'END{{print NR}}' {input.pe1})
        reads_PE2=$(bioawk -c fastx 'END{{print NR}}' {input.pe2})
        reads_PE0=0
        echo -e "{wildcards.sample}\t${{reads_PE1}}\t${{reads_PE2}}\t${{reads_PE0}}" > {output}
        """

rule summarise_count_reads:
    input:
        expand("03-data/processed_data/Zape2/{sample}.n", sample=['Zape2', 'Zape2_MinE2X', 'Zape2_PCMinE'])
    output:
        "05-results/PREP_Nextflow_EAGER_Zape2_noReads.tsv"
    message: "Summarise the number of reads per sample"
    run:
        pd.concat([pd.read_csv(fn, sep="\t", header=None, names=['sample', 'R1', 'R2', 'R0'])
                   for fn in input]) \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)

################################################################################
