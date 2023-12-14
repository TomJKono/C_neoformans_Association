#!/usr/bin/env python
"""Parse a directory of GWAS results from GRAB and/or GEMMA and compile the
results in a big table. This script will subset the results to just SNPs listed
in a file. Takes three arguments:
    1) SNP list text file, one SNP ID per line
    2) GRAB GWAS results directory
    3) GEMMA GWAS results directory"""

import os
import sys

try:
    snp_ids = sys.argv[1]
    grab_dir = sys.argv[2]
    gemma_dir = sys.argv[3]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)


def get_gwas_files(d, gwas_type):
    """Return a list of paths to GWAS output files. The 'gwas_type' argument
    specifies the default suffix to use when return a list of files. The two
    supported types are 'grab' and 'gemma.' gwas_type='grab' will return
    files ending in '.assoc' and gwas_type='gemma' will return files ending
    in '.assoc.txt'"""
    fp = os.path.abspath(os.path.expanduser(d))
    if gwas_type == 'grab':
        suffix = '.assoc'
    elif gwas_type == 'gemma':
        suffix = '.assoc.txt'
    ret_f = []
    for fname in os.listdir(fp):
        if fname.endswith(suffix):
            ret_f.append(os.path.join(fp, fname))
    return ret_f


def get_trait_name(gwas_filename, gwas_type):
    """Parse the GWAS output filename to extract the trait name. The gwas_type
    argument defines the file extension, see the get_gwas_files() function for
    which types are supported and the suffixes they correspond to."""
    if gwas_type == 'grab':
        suffix = '.assoc'
    elif gwas_type == 'gemma':
        suffix = '.assoc.txt'
    # Strip the parent path names
    bname = os.path.basename(gwas_filename)
    # Strip the suffix
    trait = bname[:-len(suffix)]
    return trait


def parse_grab(g_dir, snplist):
    """Parse through a list of GRAB output files to extract the results from
    a list of SNPs of interest."""
    grab_files = get_gwas_files(g_dir, gwas_type='grab')
    # Then, for each file, get the trait from the filename and parse the file
    # to extract results from the SNPs of interest.
    grab_dat = {}
    for gf in grab_files:
        trait = get_trait_name(gf, 'grab')
        gwas_snps = {}
        with open(gf, 'rt') as f:
            for index, line in enumerate(f):
                if index == 0:
                    continue
                tmp = line.strip().split('\t')
                snpid = tmp[0]
                # Skip SNPs that are not in the list of interest!
                if snpid not in snplist:
                    continue
                # Keep the Pvalue (6th column) and the beta value (7th column)
                gwas_snps[snpid] = (tmp[5], tmp[6])
        grab_dat[trait] = gwas_snps
    return grab_dat


def parse_gemma(g_dir, snplist):
    """Parse through a list of GEMMA output files to extract the results from
    SNPs of interest."""
    gemma_files = get_gwas_files(g_dir, gwas_type='gemma')
    gemma_dat = {}
    for gf in gemma_files:
        trait = get_trait_name(gf, 'gemma')
        gwas_snps = {}
        with open(gf, 'rt') as f:
            for index, line in enumerate(f):
                if index == 0:
                    continue
                tmp = line.strip().split('\t')
                # GEMMA reports the SNP ID as the second column
                snpid = tmp[1]
                if snpid not in snplist:
                    continue
                # Keep beta (7th column) and p_lrt (13th column), but report
                # the P-value before the beta estimate.
                gwas_snps[snpid] = (tmp[12], tmp[6])
        gemma_dat[trait] = gwas_snps
    return gemma_dat


def merge_gwas_res(g1, g2, snplist):
    """Merge two GWAS results dictionaries into a single one. We will "flip"
    the key to SNP ID, rather than trait name."""
    # Make a dictionary to hold the results
    comb_res = {s:{} for s in snplist}
    # Push the first dictionary results into the combined dict
    for trait in g1:
        for snp in snplist:
            # Get the P-value and regression beta. If we cannot find the SNP
            # in the association results, it was filtered for some reason, and
            # we will drop in NA values.
            pval, beta = g1[trait].get(snp, ('NA', 'NA'))
            comb_res[snp].update({trait: (pval, beta)})
    # Push th esecond dictionary into the combined dict
    for trait in g2:
        for snp in snplist:
            pval, beta = g2[trait].get(snp, ('NA', 'NA'))
            comb_res[snp].update({trait: (pval, beta)})
    return comb_res


def print_gwas_table(gdat):
    """Print the big matrix of GWAS results as a CSV."""
    # Identify the order in which we will print the rows (SNPs)
    snp_list = list(gdat.keys())
    snp_list = sorted(
        snp_list,
        key=lambda x: (int(x.split('_')[0]), int(x.split('_')[1])))
    # Use the first SNP to build the list of trait names for the header column
    traits = sorted(list(gdat[snp_list[0]].keys()))
    # Print the header
    header = ['SNP.ID', 'Chromosome', 'Position']
    for trt in traits:
        header.append(trt + '.Pval')
        header.append(trt + '.Beta')
    print(','.join(header))
    # Print each row
    for snp in snp_list:
        # Start with the SNP ID and then the coordinates
        coords = snp.split('_')
        toprint = [snp] + coords
        # Then append the trait Pvalues and beta values
        for trt in traits:
            toprint += list(gdat[snp][trt])
        print(','.join(toprint))
    return


def main(snps, grab, gemma):
    """Main function."""
    # Read the SNPs in as a set
    target_snps = set()
    with open(snps, 'rt') as f:
        for line in f:
            target_snps.add(line.strip())
    # Parse through the GRAB results
    grab_res = parse_grab(grab, target_snps)
    # Parse through the GEMMA results
    gemma_res = parse_gemma(gemma, target_snps)
    # Merge the two dictionaries into one big table
    comb_dat = merge_gwas_res(grab_res, gemma_res, target_snps)
    # And then print it out.
    print_gwas_table(comb_dat)
    return


main(snp_ids, grab_dir, gemma_dir)
