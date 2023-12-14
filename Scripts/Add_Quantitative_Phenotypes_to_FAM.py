#!/usr/bin/env python
"""Add the quantitative phenotypes to the PLINK FAM file for analysis with
GEMMA. We will hard-code the phenotype column index values for the sake of
convenience of development. Takes two arguments:
    1) PLINK FAM file
    2) Phenotypes CSV file"""

import sys

try:
    fam_in = sys.argv[1]
    pheno_in = sys.argv[2]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)


# Define the column offsets (in Python counting, from 0) for the quantitative
# phenotypes for GEMMA. These are the following columns (in order):
#   Virulence..dpi.
#   IL.4
#   IL.5
#   IL.13
#   IFN.γ
#   IL.1b
#   IL.2
#   IL.6
#   IL.12p70
#   GM.CSF
#   TNF.α
#   IL.18
#   Brain.CFU
#   Lung.CFU
PHENO_COLS = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 30, 31]


def parse_phenotype(pheno_csv):
    """Parse the phenotypic CSV file and return a dictionary with the
    quantitative phenotypes in it. The key will be the within-family ID, and
    the value will be a list of phenotypic values."""
    p_dict = {}
    with open(pheno_csv, 'rt') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            tmp = line.strip().split(',')
            isolate = tmp[0]
            i_phenos = [tmp[i] for i in PHENO_COLS]
            p_dict[isolate] = i_phenos
    return p_dict


def main(fam, pheno):
    """Main function."""
    pheno_dat = parse_phenotype(pheno)
    with open(fam, 'rt') as f:
        for line in f:
            tmp = line.strip().split()
            iid = tmp[1]
            i_phen = pheno_dat.get(iid)
            toprint = tmp[:-1] + i_phen
            print(' '.join(toprint))
    return


main(fam_in, pheno_in)
