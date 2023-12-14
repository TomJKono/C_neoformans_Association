#!/bin/bash
# Run the GEMMA analysis of the quantitative traits.

# Path to GEMMA executable
GEMMA="/storage/Data/Sci-Software/gemma-0.98.5/gemma-0.98.5-linux-static-AMD64"

# Path to the PLINK file basename
#   This is for the full dataset
# BFILE_BASE="/home/tom/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Data/38_Strain_PLINK/MAF_Filtered/Cn_38_Strains.MAF"
#   This is for the ST93A subset
BFILE_BASE="/home/tom/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Data/38_Strain_PLINK/ST93A_Isolates/ST93A_Isolates"


# Output directory
#   For the full dataset
# OUT_DIR="/home/tom/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/Associations"
OUT_DIR="/home/tom/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/ST93A_Isolates/Associations/GEMMA"
cd "${OUT_DIR}"

# Generate the GRM:
#   For the full dataset
# "${GEMMA}" -bfile "${BFILE_BASE}" -gk 1 -o "Cn_38_Strains.GRM"
#   For the ST93A subset
"${GEMMA}" -bfile "${BFILE_BASE}" -gk 1 -o "ST93A.GRM"

# Based on the GRM basename, get the full path to the GRM file
#   For the full dataset
# GRM="/home/tom/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/Associations/output/Cn_38_Strains.GRM.cXX.txt"
#   For the ST93A subset
GRM="/home/tom/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/ST93A_Isolates/Associations/GEMMA/output/ST93A.GRM.cXX.txt"

# Then run the associations in a loop. This will be kind of ugly...
#   First, define an array of phenotype names for output files
PHENOTYPES=(
    "Virulence_dpi"
    "IL-4"
    "IL-5"
    "IL-13"
    "IFN-Gamma"
    "IL-1b"
    "IL-2"
    "IL-6"
    "IL-12p70"
    "GM-CSF"
    "TNF-Alpha"
    "IL-18"
    "Brain_CFU"
    "Lung_CFU")

for i in {0..13}
do
    # Get the current phenotype name from the phenotype array
    CURR_PHENO="${PHENOTYPES[${i}]}"
    # Calculate the offset for the phenotype in the FAM file (i+1)
    PHENO_IDX=$((${i} + 1))
    # Run the GEMMA analysis
    "${GEMMA}" -bfile "${BFILE_BASE}" -n "${PHENO_IDX}" -k "${GRM}" -o "${CURR_PHENO}" -miss 0.50 -maf 0.02 -r2 1.0 -hwe 0 -lmm 4
done
