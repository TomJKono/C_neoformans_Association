#!/bin/bash
# Simple script to isolate SNPs in genes of interest.

# Define a regular expression to 'grep' out the genes of interest
GENE_RE=$'(CNAG_07528|CNAG_02112|CNAG_05662)'
# Define the input file
FILE_IN="Data/38_Strain_Orig/652 Varients.csv"
# Define the output file
FILE_OUT="Resources/38_Strain_SNPs_in_GOI.txt"

grep -E "${GENE_RE}" "${FILE_IN}" | awk -F ',' '{print $1"_"$2}' > "${FILE_OUT}"
