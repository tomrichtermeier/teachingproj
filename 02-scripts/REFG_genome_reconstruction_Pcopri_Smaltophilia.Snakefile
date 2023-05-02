################################################################################
# Project: Training project
# Part: Reference-based genome reconstruction
# Step: Reconstruction of P. copri and S. maltophilia
#
# Dependent on:
#   - QUAL_damageprofiles.Snakefile
#
# Alex Huebner, 27/04/23
################################################################################

import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("03-data/processed_data/{sample}_1.fastq.gz")
################################################################################

wildcard_constraints:
    sample = "[ES]RS[0-9]+"

rule all:
    input:
        expand("04-analysis/refgenome_reconst/{sample}.freebayes_loose.fasta.gz", sample=SAMPLES),
        expand("04-analysis/refgenome_reconst/{sample}.freebayes_conserv.fasta", sample=[s for s in SAMPLES if s not in ['ERS11453815', 'ERS11453816']]),
        expand("04-analysis/refgenome_reconst/{sample}.majority.vcf.gz", sample=SAMPLES)

#### Prepare reference for genotyping with freeBayes ###########################

rule uncompress_reffasta:
    output:
        temp("tmp/genome_reconst/{genome}.fasta")
    message: "Decompress the FastA file of the reference genome: {wildcards.genome}"
    resources:
        mem = 4,
        cores = 1
    params:
        fasta = "tmp/damageprofiles/{genome}.fna.gz"
    shell:
        """
        gunzip -c {params.fasta} > {output}
        """

rule faidx_reffasta:
    input:
        "tmp/genome_reconst/{genome}.fasta"
    output:
        temp("tmp/genome_reconst/{genome}.fasta.fai")
    message: "Generate FastA index for the reference genome: {wildcards.genome}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 4,
        cores = 1
    shell:
        """
        samtools faidx {input}
        """

################################################################################

#### Genotyping using freeBayes ################################################

rule freebayes:
    input:
        fa = lambda wildcards: "tmp/genome_reconst/Pcopri.fasta" if wildcards.sample[:3] == "SRS" or wildcards.sample[:6] == "ERS418" else "tmp/genome_reconst/Smaltophilia.fasta",
        fai = lambda wildcards: "tmp/genome_reconst/Pcopri.fasta.fai" if wildcards.sample[:3] == "SRS" or wildcards.sample[:6] == "ERS418" else "tmp/genome_reconst/Smaltophilia.fasta.fai",
    output:
        pipe("tmp/genome_reconst/{sample}.vcf")
    message: "Genotype the contigs using freeBayes in parallel mode: {wildcards.sample}"
    conda: "ENVS_freebayes.yaml"
    group: "freebayes"
    resources:
        mem = 8,
        cores = 1
    params:
        bam = "04-analysis/damageprofiles/{sample}.calmd.bam"
    threads: 1
    shell:
        """
        freebayes -f {input.fa} \
            --report-monomorphic \
            -C 1 -F 0.05 -p 1 \
            --haplotype-length 0 \
            -q 30 -m 20 {params.bam} > {output}
        """

rule compress_vcf:
    input:
        "tmp/genome_reconst/{sample}.vcf"
    output:
        "04-analysis/refgenome_reconst/{sample}.freebayes.vcf.gz"
    message: "Compress the VCF file produced by freebayes: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    group: "freebayes"
    resources:
        mem = 4,
        cores = 1
    threads: 1
    shell:
        """
        bgzip -c {input} > {output}
        """

rule bcftools_filter:
    input:
        "04-analysis/refgenome_reconst/{sample}.freebayes.vcf.gz"
    output:
        vcf = "04-analysis/refgenome_reconst/{sample}.filter.vcf.gz",
        tbi = "04-analysis/refgenome_reconst/{sample}.filter.vcf.gz.tbi"
    message: "Discard low-quality differences between MEGAHIT and freebayes consensus: {wildcards.sample}"
    conda: "ENVS_bcftools.yaml"
    resources:
        mem = 4,
        cores = 1
    shell:
        """
        bcftools view \
            -v snps,mnps \
            -i 'QUAL >= 30 || (QUAL >= 20 && INFO/AO >= 3)' {input} | \
        bgzip > {output.vcf}
        bcftools index -t {output.vcf}
        """

################################################################################

#### Majority calling ##########################################################

rule download_bamcaller:
    output:
        "02-scripts/pyscripts/bam-caller.py"
    message: "Download script for bam-caller"
    params:
        url = "https://raw.githubusercontent.com/bodkan/bam-caller/master/bam-caller.py"
    shell:
        "wget -O {output} {params.url} && chmod u+x {output}"

rule bamcaller:
    input:
        "02-scripts/pyscripts/bam-caller.py"
    output:
        temp("04-analysis/refgenome_reconst/{sample}.majority.vcf")
    message: "Infer majority alleles with bam-caller: {wildcards.sample}"
    conda: "ENVS_bamcaller.yaml"
    group: "bamcaller"
    params:
        bam = "04-analysis/damageprofiles/{sample}.calmd.bam",
        prefix = "04-analysis/refgenome_reconst/{sample}.majority"
    shell:
        """
        {input} --bam {params.bam} \
            --strategy majority \
            --proportion 0.67 \
            --mincov 1 \
            --minbq 30 \
            --minmq 20 \
            --sample-name {wildcards.sample} \
            --output {params.prefix}
        """

rule compress_bamcaller:
    input:
        "04-analysis/refgenome_reconst/{sample}.majority.vcf"
    output:
        "04-analysis/refgenome_reconst/{sample}.majority.vcf.gz"
    message: "Compress the VCF file produced by bamcaller: {wildcards.sample}"
    conda: "ENVS_samtools.yaml"
    group: "bamcaller"
    resources:
        mem = 4,
        cores = 1
    threads: 1
    shell:
        """
        bgzip -c {input} > {output}
        """

################################################################################

#### Call consensus ############################################################

rule bcftools_consensus:
    input:
        vcf = "04-analysis/refgenome_reconst/{sample}.filter.vcf.gz",
        tbi = "04-analysis/refgenome_reconst/{sample}.filter.vcf.gz.tbi"
    output:
        "04-analysis/refgenome_reconst/{sample}.freebayes_loose.fasta.gz"
    message: "Correct the consensus sequence of the contigs: {wildcards.sample}"
    conda: "ENVS_bcftools.yaml"
    resources:
        mem = 8,
        cores = 2
    params:
        fasta = lambda wildcards: "tmp/damageprofiles/Pcopri.fna.gz" if wildcards.sample[:3] == "SRS" or wildcards.sample[:6] == "ERS418" else "tmp/damageprofiles/Smaltophilia.fna.gz"
    threads: 2
    shell:
        """
        cat {params.fasta} | bcftools consensus {input.vcf} | bgzip > {output}
        """

rule vcf2fasta_freebayes:
    input:
        "04-analysis/refgenome_reconst/{sample}.freebayes.vcf.gz"
    output:
        "04-analysis/refgenome_reconst/{sample}.freebayes_conserv.fasta"
    message: "Convert the freeBayes VCF file into FastA file: {wildcards.sample}"
    conda: "ENVS_vcf2fasta.yaml"
    resources:
        mem = 8,
        cores = 1
    params:
        fasta = lambda wildcards: "tmp/damageprofiles/Pcopri.fna.gz" if wildcards.sample[:3] == "SRS" or wildcards.sample[:6] == "ERS418" else "tmp/damageprofiles/Smaltophilia.fna.gz"
    shell:
        """
        02-scripts/pyscripts/vcf2fasta.py \
            -i {input} \
            -o {output} \
            -r {params.fasta} \
            --minqual_fallback 20 --mincov_fallback 3
        """

################################################################################
