#!/usr/bin/env python
"""Generate SNP IDs for the SNPs of interest using a gene-chromosome table and
the coordinates from Katrina's table. Takes two arguments:
    1) gene-chromosome table
    2) SNPs of interest CSV"""

import sys

try:
    gene_chr = sys.argv[1]
    snp_csv = sys.argv[2]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)


def main(gc_table, snp_table):
    """Main function."""
    # First, parse the gene->chromosome file and store it as a dict
    gene_chr_dict = {}
    with open(gc_table, 'rt') as f:
        for line in f:
            tmp = line.strip().split(',')
            # First column is chromosome, second is gene name
            gene_chr_dict[tmp[1]] = tmp[0]
    # Then iterate through the SNP table and build the SNP ID
    with open(snp_table, 'rt') as f:
        for index, line in enumerate(f):
            if index == 0:
                header = line.strip().split(',')
                header = ['SNP_ID'] + header
                print(','.join(header))
            else:
                tmp = line.strip().split(',')
                coord = tmp[0]
                gene = tmp[1]
                chrom = gene_chr_dict.get(gene)
                snp_id = chrom + '_' + coord
                toprint = [snp_id] + tmp
                print(','.join(toprint))
    return


main(gene_chr, snp_csv)
