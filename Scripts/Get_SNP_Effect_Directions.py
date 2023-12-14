#!/usr/bin/env python
"""Take a table of SNPs and traits of interest and report the directionality
of the effect substition (trait up/down) for each trait for each SNP. Takes
three arguments:
    1) SNP-of-interest and trait table
    2) GEMMA output dir
    3) GRAB output dir"""

import sys
import os
import pprint

try:
    snp_table = sys.argv[1]
    gemma_dir = sys.argv[2]
    grab_dir = sys.argv[3]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)


def parse_grab(d):
    """Parse a directory of GRAB output files and return a nested dictionary
    where the top-level key is the trait name and the value is a dictionary of
    SNP IDs and effect directions."""
    grab_fp = os.path.abspath(os.path.expanduser(d))
    grab_files = [
        os.path.join(grab_fp, f)
        for f
        in os.listdir(grab_fp)
        if f.endswith('.assoc')]
    # Then go through this list of files and build up the dictionary
    grab_results = {}
    for grab_fname in grab_files:
        bname = os.path.basename(grab_fname)
        # Strip off the .assoc and replace underscores with spaces to get the
        # human-readable trait name
        trait = bname.replace('.assoc', '').replace('_', ' ')
        trait_res = {}
        with open(grab_fname, 'rt') as f:
            for index, line in enumerate(f):
                if index == 0:
                    continue
                tmp = line.strip().split('\t')
                snp_id = tmp[0]
                pval = tmp[5]
                beta = tmp[6]
                # set the direction based on the sign of the beta
                try:
                    beta = float(beta)
                    if beta > 0:
                        direction = 'Positive'
                    elif beta < 0:
                        direction = 'Negative'
                    elif beta == 0:
                        direction = 'No Effect'
                except TypeError:
                    direction = 'Missing'
                trait_res[snp_id] = (str(beta), pval, direction)
        grab_results[trait] = trait_res
    return grab_results


def parse_gemma(d):
    """Parse a directory of GEMMA output files and return a nested dictionary
    where the top-level key is the trait name and the lower key is the SNP ID.
    """
    gemma_fp = os.path.abspath(os.path.expanduser(d))
    gemma_files = [
        os.path.join(gemma_fp, f)
        for f
        in os.listdir(gemma_fp)
        if f.endswith('.assoc.txt')]
    # Build up the dictionary in the same way as for the GRAB files
    gemma_results = {}
    for gemma_fname in gemma_files:
        bname = os.path.basename(gemma_fname)
        # Strip off the .assoc and replace underscores with spaces to get the
        # human-readable trait name
        trait = bname.replace('.assoc.txt', '').replace('_', ' ')
        trait_res = {}
        with open(gemma_fname, 'rt') as f:
            for index, line in enumerate(f):
                if index == 0:
                    continue
                tmp = line.strip().split('\t')
                snp_id = tmp[1]
                pval = tmp[13]
                beta = tmp[7]
                # set the direction based on the sign of the beta
                try:
                    beta = float(beta)
                    if beta > 0:
                        direction = 'Positive'
                    elif beta < 0:
                        direction = 'Negative'
                    elif beta == 0:
                        direction = 'No Effect'
                except TypeError:
                    direction = 'Missing'
                trait_res[snp_id] = (str(beta), pval, direction)
        gemma_results[trait] = trait_res
    return gemma_results


def main(snps, gemma_d, grab_d):
    """Main function."""
    # First, parse out the GRAB results
    grab_dat = parse_grab(grab_d)
    # Do the same for the GEMMA results
    gemma_dat = parse_gemma(gemma_d)
    # OK, now we can iterate through the SNP table and get the directions.
    # Print a header
    print('SNP.ID,Gene,Coord,Trait,Beta,P.Value,Direction')
    with open(snps, 'rt') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            tmp = line.strip().split(',')
            snp_id = tmp[0]
            coord = tmp[1]
            gene = tmp[2]
            traits = tmp[3:]
            # Remove empty entries from the traits list
            traits = [x for x in traits if x != '']
            for tr in traits:
                # Replace CFUs with CFU
                tr = tr.replace('CFUs', 'CFU')
                # Replace 'IFNg' with 'IFN-Gamma'
                tr = tr.replace('IFNg', 'IFN-Gamma')
                # Replace 'TNF alpha' with 'TNF-Alpha'
                tr = tr.replace('TNF alpha', 'TNF-Alpha')
                # Replace '37C' with '37'
                tr = tr.replace('37C', '37')
                # Replace 'L-DOPA 30C' with 'L-DOPA-30'
                tr = tr.replace('L-DOPA 30C', 'L-DOPA-30')
                # Replace 'Virlence category' with 'Virulence Category'
                tr = tr.replace('Virulence category', 'Virulence Category')
                if tr in grab_dat:
                    beta, pval, direction = grab_dat[tr].get(
                        snp_id, ('NA', 'NA', 'Missing'))
                elif tr in gemma_dat:
                    beta, pval, direction = gemma_dat[tr].get(
                        snp_id, ('NA', 'NA', 'Missing'))
                else:
                    beta, direction = ('', '', 'Cannot find trait')
                toprint = [snp_id, gene, coord, tr, beta, pval, direction]
                print(','.join(toprint))
    return


main(snp_table, gemma_dir, grab_dir)
