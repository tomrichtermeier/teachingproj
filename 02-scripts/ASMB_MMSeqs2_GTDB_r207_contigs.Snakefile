################################################################################
# Project: Training project
# Part: De novo assembly
# Step: Assign all contigs longer than 500 bp to a GTDB taxonomy using MMSeqs2
# 
# Alex Huebner, 12/04/23
################################################################################

from glob import glob
import gzip
import os

import pandas as pd
import pyfastx

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("04-analysis/metagenome_assembly/alignment/metaspades/{sample}-metaspades.fasta.gz")
################################################################################

#### SPLITS FOR CONTIG SIZES ###################################################
SPLITS = {0 : (500, 1000),
          1 : (1000, 2000),
          2 : (2000, 5000),
          3 : (5000, 10000),
          4 : (10000, 25000),
          5 : (25000, 50000),
          6 : (50000, 5000000)}
################################################################################

wildcard_constraints:
    split = "[0-9]"

rule all:
    input:
         "05-results/ASMB_MMSeqs2_GTDB_taxassignment_contigs.tsv.gz"

rule concat_contigs:
    output:
        temp("tmp/mmseqs2/{split}.fasta.gz")
    message: "Concatenate all contigs of split {wildcards.split} into a single FastA"
    resources:
        mem = 4,
        cores = 1
    run:
        low = SPLITS[int(wildcards.split)][0]
        high = SPLITS[int(wildcards.split)][1]
        with gzip.open(output[0], "wt") as outfile:
            for sample in SAMPLES:
                for name, seq in pyfastx.Fasta(f"04-analysis/metagenome_assembly/alignment/metaspades/{sample}-metaspades.fasta.gz", build_index=False):
                    if len(seq) >= low and len(seq) < high:
                        outfile.write(f">{sample}:{name}\n{seq}\n")

rule create_db:
    input:
        "tmp/mmseqs2/{split}.fasta.gz"
    output:
        temp("tmp/mmseqs2/{split}.contigs")
    message: "Create database of contigs: {wildcards.split}"
    resources:
        mem = 16,
        cores = 1
    wrapper:
        "file:///home/alexander_huebner/github/snakemake-wrappers/bio/mmseqs2/createdb"

rule install_refdb:
    output:
        "03-data/refdbs/mmseqs2_gtdb_r207_db/mmseqs2_gtdb_r207_db"
    message: "Download and install the MMSeqs2 GTDB r207"
    resources:
        mem = 24,
        cores = 1
    params:
        dbname = "GTDB",
        prefix = "03-data/refdbs/mmseqs2_gtdb_r207_db/mmseqs2_gtdb_r207_db"
    wrapper:
        "file:///home/alexander_huebner/github/snakemake-wrappers/bio/mmseqs2/download_database"

rule taxonomy:
    input:
        contigs = "tmp/mmseqs2/{split}.contigs",
        db = "03-data/refdbs/mmseqs2_gtdb_r207_db/mmseqs2_gtdb_r207_db"
    output:
        "tmp/mmseqs2/{split}.mmseqs2_gtdb.index"
    message: "Assign taxonomy via the GTDB for contigs: {wildcards.split}"
    resources:
        mem = 750,
        cores = 32
    params:
        prefix = "tmp/mmseqs2/{split}.mmseqs2_gtdb",
        extra = "--search-type 2 -s 5.0 --orf-filter-s 4.0 --lca-ranks kingdom,phylum,class,order,family,genus,species --tax-lineage 1"
    threads: 32
    wrapper:
        "file:///home/alexander_huebner/github/snakemake-wrappers/bio/mmseqs2/taxonomy"

rule create_tsv:
    input:
        contigs = "tmp/mmseqs2/{split}.contigs",
        assignments = "tmp/mmseqs2/{split}.mmseqs2_gtdb.index"
    output:
        temp("tmp/mmseqs2/{split}.mmseqs2_gtdb.tsv")
    message: "Convert MMSeqs2 GTDB results to TSV: {wildcards.split}"
    resources:
        mem = 48,
        cores = 1
    params:
        query = "tmp/mmseqs2/{split}.contigs",
        target = "",
        prefix = "tmp/mmseqs2/{split}.mmseqs2_gtdb"
    wrapper:
        "file:///home/alexander_huebner/github/snakemake-wrappers/bio/mmseqs2/createtsv"

rule mmseqs2_annotatetsv:
    input:
        "tmp/mmseqs2/{split}.mmseqs2_gtdb.tsv"
    output:
        "tmp/mmseqs2/{split}.mmseqs2_gtdb.annot.tsv"
    message: "Add header to MMSeqs2 table: {wildcards.split}"
    run:
        pd.read_csv(input[0], sep="\t", header=None,
                    usecols=list(range(8)) + [9],
                    names=['contig', 'NCBItaxID', 'NCBIrank', 'NCBItaxName',
                        'nFrags', 'retainedFrags', 'taxassignedFrags',
                        'fractionAgreement', 'lineage']) \
            .to_csv(output[0], sep="\t", index=False)

rule combine_mmseqs2_tsvs:
    input:
        expand("tmp/mmseqs2/{split}.mmseqs2_gtdb.annot.tsv", split=SPLITS)
    output:
        "05-results/ASMB_MMSeqs2_GTDB_taxassignment_contigs.tsv.gz"
    message: "Concatenate the MMSeqs2 taxonomy results"
    run:
        tbl = pd.concat([pd.read_csv(fn, sep="\t")
                   for fn in input])
        tbl['sample'] = tbl['contig'].str.split(":").str[0]
        tbl['contigID'] = tbl['contig'].str.extract(r'[A-Z]+[0-9]+:NODE_([0-9]+)_length_[0-9]+_cov_[0-9\.]+').astype(int)
        tbl['contig'] = tbl['contig'].str.split(":").str[1]

        tbl.sort_values(['sample', 'contigID']) \
            .drop(['contigID'], axis=1) \
            .to_csv(output[0], sep="\t", index=False, compression="gzip", float_format="%.2f")
