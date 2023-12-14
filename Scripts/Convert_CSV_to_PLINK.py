#!/usr/bin/env python
"""Convert the genotype matrix CSV to PLINK PED/MAP file set. Takes two
arguments:
    1) CSV export of Excel genotype matrix
    2) Base name for the PLINK file set"""

import sys

try:
    csv_in = sys.argv[1]
    bname_out = sys.argv[2]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)

marker_info = []
sample_info = {}

with open(csv_in, 'rt') as f:
    for index, line in enumerate(f):
        if index == 0:
            header = line.strip().split(',')
            snames = header[9:]
            # Initialize the dictionary for genotype info
            for sn in snames:
                # Skip the "CP" column - not sure what this is
                if sn == 'CP':
                    continue
                sample_info[sn] = []
        else:
            tmp = line.strip().split(',')
            chrom = tmp[0]
            bp = tmp[1]
            # If there are indels, then we will treat them a little
            # differently: indels will be recoded to A/G, with 'A' denoting
            # the short alelle and 'G' denoting the long allele.
            ref = tmp[2]
            alt = tmp[3]
            indel = bool(len(ref) > 1 or len(alt) > 1)
            # make up a SNP name based on the position
            snp_name = chrom + '_' + bp
            # store the marker info in the list that will be used to
            # build the MAP file. The MAP file fields are:
            #   1: Chromosome
            #   2: SNP identifier
            #   3: Genetic position in centimorgans ('0' for unknown)
            #   4: Physical position
            marker_info.append((chrom, snp_name, '0', bp))
            # Add the genotypes into the dictionary that will be used to build
            # the PED file
            genos = tmp[9:]
            # Get the alleles of the marker
            alleles = set(genos)
            for sn, gt in zip(snames, genos):
                if sn == 'CP':
                    continue
                # Set the missing code properly
                if gt == 'NA':
                    new_gt = '0'
                else:
                    # If we are dealing with an indel, then we have to check
                    # the length of the call
                    if indel:
                        if len(gt) == 1:
                            new_gt = 'A'
                        else:
                            new_gt = 'G'
                    else:
                        new_gt = gt
                sample_info[sn].append(new_gt)

# Write the MAP file
map_out = open(bname_out + '.map', 'wt')
for row in marker_info:
    map_out.write('\t'.join(row) + '\n')
map_out.flush()
map_out.close()

# Write the PED file:
ped_out = open(bname_out + '.ped', 'wt')
for samplename in sorted(sample_info):
    gt_row = []
    fid = '0'
    iid = samplename
    pid = '0'
    mid = '0'
    sex = '0'
    pheno = '-9'
    gt_row = [fid, iid, pid, mid, sex, pheno]
    for gt in sample_info[samplename]:
        # We need two append statements here because even though the fungi
        # are haploid, PLINK files are written for diploid data
        gt_row.append(gt)
        gt_row.append(gt)
    # print out the row
    ped_out.write('\t'.join(gt_row) + '\n')
ped_out.flush()
ped_out.close()

