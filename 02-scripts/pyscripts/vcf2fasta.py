#!/usr/bin/env python
########################################################################
# Convert VCF format to FastA without assuming REF as default allele
#
# Alex Huebner, 27/04/23
########################################################################

from __future__ import print_function
import argparse
import sys

import allel
import numpy as np
import pyfastx


def main():
    ''' Convert VCF file into FastA sequence
    '''
    # Read VCF data
    vcf = allel.read_vcf(Args['input'],
                         fields=['variants/CHROM', 'variants/POS',
                                 'variants/REF', 'variants/ALT',
                                 'variants/QUAL', 'calldata/*'],
                         numbers={'calldata/GT': Args['ploidy']})

    if vcf is None:
        print("This VCF file contains no entries.", file=sys.stderr)
        sys.exit(1)

    # Create default sequences with Ns for each contig
    contigs = {contig: len(seq)
               for contig, seq in pyfastx.Fasta(Args['ref'], build_index=False)}
    consensus = {contig: ['N'] * contigs[contig]
                 for contig in contigs}

    # Ref genotype
    if Args['ploidy'] == 1:
        refgt = 0

    # Replace Ns with genotypes when data fulfills quality requirements
    for contig in contigs:
        if contig in vcf['variants/CHROM']:
            sites_idx = np.where(vcf['variants/CHROM'] == contig)[0]
            for i in sites_idx:
                if ((vcf['variants/QUAL'][i] >= Args['minqual']
                     and vcf['calldata/DP'][i][0] > Args['mincov'])
                    or (vcf['variants/QUAL'][i] >= Args['minqual_fallback']
                        and vcf['calldata/DP'][i][0] > Args['mincov_fallback'])
                    or (vcf['calldata/DP'][i][0] > Args['mincov'] and
                        vcf['calldata/GT'][i][0] == refgt)):  # reference allele
                    if vcf['calldata/GT'][i][0] == 0:
                        gt = vcf['variants/REF'][i]
                    elif vcf['calldata/GT'][i][0] == 1:
                        gt = vcf['variants/ALT'][i][0]
                    else:
                        continue
                    consensus[contig][vcf['variants/POS'][i] - 1] = gt

    # Write to FastA
    with open(Args['output'], "wt") as outfile:
        for contig in contigs.keys():
            outfile.write(f">{contig}\n")
            outfile.write("".join(consensus[contig]) + "\n")


# Argument parser
Parser = argparse.ArgumentParser(description='Converts VCF to FastA without '
                                 'assuming the REF allele as default')
Parser.add_argument('-i', '--input', required=True, help='VCF file')
Parser.add_argument('-o', '--output', required=True, help='FastA file')
Parser.add_argument('-r', '--ref', required=True, help='reference FastA file')
Parser.add_argument('--ploidy', default=1, type=int,
                    help='ploidy used for genotyping [1]')
Parser.add_argument('--minqual', default=30, type=int,
                    help='minimal genotype quality [30]')
Parser.add_argument('--mincov', default=1, type=int,
                    help='minimal coverage [1]')
Parser.add_argument('--minqual_fallback', default=30, type=int,
                    help='minimal genotype quality to consider when more '
                    'coverage is available [30]')
Parser.add_argument('--mincov_fallback', default=1, type=int,
                    help='minimal coverage required for fallback evaluation [1]')
Args = vars(Parser.parse_args())

if __name__ == '__main__':
    main()
