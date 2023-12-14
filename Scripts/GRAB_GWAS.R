# Use the GRAB package to run a GWAS on each of the traits, individually.

library(GRAB)

# Paths to input files
gh_base <- "/Users/tomkono/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Data/38_Strain_PLINK/"
phenotypes <- paste(gh_base, "All_Phenotypes.csv", sep="")
#   This is the full dataset
# genotypes <- paste(gh_base, "MAF_Filtered/Cn_38_Strains.MAF.bed", sep="")
#   This is the ST93A subset
genotypes <- paste(gh_base, "ST93A_Isolates/ST93A_Isolates.bed", sep="")
#   For the full dataset
# pcs <- paste(gh_base, "MAF_Filtered/Cn_38_Strains.MAF.PCA.eigenvec", sep="")
#   For the ST93A subset
pcs <- paste(gh_base, "ST93A_Isolates/ST93A_Isolates.PCA.eigenvec", sep="")

# Paths to output file
# out_base <- "/Users/tomkono/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/Associations/"
out_base <- "/Users/tomkono/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/ST93A_Isolates/Associations/"

# Read in the data files
phen_dat <- read.csv(phenotypes, header=TRUE)
pc_dat <- read.table(pcs, header=FALSE)

# Set the column names of the phenotypes to easy-to-work-with names.
colnames(phen_dat) <- c("Isolate", "Vir.DPI", "Vir.Cat", "IL.4", "IL.5", "IL.13",
                       "IFN.Gamma", "IL.1b", "IL.2", "IL.6", "IL.12p70", "GM.CSF",
                       "TNF.Alpha", "IL.18", "X37", "pH.8.15", "NaCl", "Congo.Red",
                       "SDS", "Caffeine", "CFW", "Niger.Seed.25", "Niger.Seed.30",
                       "L.DOPA.25", "L.DOPA.30", "Urease", "Total.Cell",
                       "Cell.Body", "Capsule", "Titan", "Brain.CFU", "Lung.CFU")

# Subset the phenotype data to just the samples we are processing
sample_list <- read.table("/Users/tomkono/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Resources/ST93A_Isolates.txt", header=FALSE)$V2
phen_dat <- phen_dat[phen_dat$Isolate %in% sample_list,]

# Append PC1 and PC2
phen_dat$PC1 <- pc_dat$V3
phen_dat$PC2 <- pc_dat$V4

# Cast the appropriate ordinal categorical variables to factor and set the
# levels. These are the ordinal categorical variables:
#   Virulence category
#   Growth rate with "37"
#   Growth rate with NaCl
#   Growth rate with Congo Red
#   Growth rate with SDS
#   Growth rate with caffeine
#   Growth rate wtih CFW
#   Growth rate with Niger Seed 25
#   Growth rate with Niger Seed 30
#   Growth rate with L-DOPA 25
#   Growth rate with L-DOPA 30
#   Growth rate wtih Urease
#   Total Cell
#   Cell Body
#   Capsule
#   Titan
phen_dat$Vir.Cat[phen_dat$Vir.Cat == "Latent"] <- 0
phen_dat$Vir.Cat[phen_dat$Vir.Cat == "Typical non-CNS"] <- 1
phen_dat$Vir.Cat[phen_dat$Vir.Cat == "Typical CNS"] <- 2
phen_dat$Vir.Cat[phen_dat$Vir.Cat == "Hypervirulent"] <- 3
phen_dat$Vir.Cat <- factor(phen_dat$Vir.Cat, levels=c(0, 1, 2, 3))
phen_dat$X37 <- factor(phen_dat$X37, levels=c(-2, -1, 0, 1))
phen_dat$pH.8.15 <- factor(phen_dat$pH.8.15, levels=c(-3, -2, -1, 0))
phen_dat$NaCl <- factor(phen_dat$NaCl, levels=c(-3, -2, -1, 0, 1))
phen_dat$Congo.Red <- factor(phen_dat$Congo.Red, levels=c(-5, -3, -2, -1, 0, 1, 3))
phen_dat$SDS <- factor(phen_dat$SDS, levels=c(-5, -4, -3, -2))
phen_dat$Caffeine <- factor(phen_dat$Caffeine, levels=c(-4, -3, -2, -1, 0, 1, 2))
phen_dat$CFW <- factor(phen_dat$CFW, levels=c(-5, -3, -2, -1, 0, 1, 2))
phen_dat$Niger.Seed.25 <- factor(phen_dat$Niger.Seed.25, levels=c(1, 2, 3, 4, 5))
phen_dat$Niger.Seed.30 <- factor(phen_dat$Niger.Seed.30, levels=c(0, 1, 2, 3, 4, 5))
phen_dat$L.DOPA.25 <- factor(phen_dat$L.DOPA.25, levels=c(1, 2, 3, 4, 5))
phen_dat$L.DOPA.30 <- factor(phen_dat$L.DOPA.30, levels=c(1, 2, 3, 4))
phen_dat$Urease <- factor(phen_dat$Urease, levels=c(-1, 0, 1))
phen_dat$Total.Cell <- factor(phen_dat$Total.Cell, levels=c(-4, -3, -2, -1, 2))
phen_dat$Cell.Body <- factor(phen_dat$Cell.Body, levels=c(-1, 1, 2))
phen_dat$Capsule <- factor(phen_dat$Capsule, levels=c(-3, -2, -1, 0))
phen_dat$Titan <- factor(phen_dat$Titan, levels=c(0, 1))
# The others are true quantitative traits, which can be analyzed with GEMMA

# This is going to be a big ugly block of repetitive statements!! 17
# phenotypes. Yikes!
gwas_control <- list(
    showInfo=FALSE,
    LOCO=FALSE,
    tolTau=1e-3,
    tolBeta=1e-3,
    maxiter=10000,
    min_mac_marker=0,
    min_maf_marker=0,
    missing_cutoff=0.5)

vir_cat.res <- GRAB.NullModel(
    formula=Vir.Cat~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
x37.res <- GRAB.NullModel(
    formula=X37~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
pH.8.15.res <- GRAB.NullModel(
    formula=pH.8.15~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
nacl.res <- GRAB.NullModel(
    formula=NaCl~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
congo_red.res <- GRAB.NullModel(
    formula=Congo.Red~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
sds.res <- GRAB.NullModel(
    formula=SDS~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
caffeine.res <- GRAB.NullModel(
    formula=Caffeine~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
cfw.res <- GRAB.NullModel(
    formula=CFW~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
niger_seed_25.res <- GRAB.NullModel(
    formula=Niger.Seed.25~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
niger_seed_30.res <- GRAB.NullModel(
    formula=Niger.Seed.30~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
l_dopa_25.res <- GRAB.NullModel(
    formula=L.DOPA.25~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
l_dopa_30.res <- GRAB.NullModel(
    formula=L.DOPA.30~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
urease.res <- GRAB.NullModel(
    formula=Urease~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
total_cell.res <- GRAB.NullModel(
    formula=Total.Cell~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
cell_body.res <- GRAB.NullModel(
    formula=Cell.Body~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
capsule.res <- GRAB.NullModel(
    formula=Capsule~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)
titan.res <- GRAB.NullModel(
    formula=Titan~PC1+PC2, data=phen_dat, subjData=phen_dat$Isolate,
    method="POLMM", traitType="ordinal", GenoFile=genotypes, control=gwas_control)

# Then print the marker-level results for these association tests
marker_control <- list(min_mac_marker=0, min_maf_marker=0, missing_cutoff=0.5)

GRAB.Marker(vir_cat.res, GenoFile=genotypes, OutputFile=paste(out_base, "Virulence_Category.assoc", sep=""), control=marker_control)
GRAB.Marker(x37.res, GenoFile=genotypes, OutputFile=paste(out_base, "37.assoc", sep=""), control=marker_control)
GRAB.Marker(pH.8.15.res, GenoFile=genotypes, OutputFile=paste(out_base, "pH-8.15.assoc", sep=""), control=marker_control)
GRAB.Marker(nacl.res, GenoFile=genotypes, OutputFile=paste(out_base, "NaCl.assoc", sep=""), control=marker_control)
GRAB.Marker(congo_red.res, GenoFile=genotypes, OutputFile=paste(out_base, "Congo_Red.assoc", sep=""), control=marker_control)
GRAB.Marker(sds.res, GenoFile=genotypes, OutputFile=paste(out_base, "SDS.assoc", sep=""), control=marker_control)
GRAB.Marker(caffeine.res, GenoFile=genotypes, OutputFile=paste(out_base, "Caffeine.assoc", sep=""), control=marker_control)
GRAB.Marker(cfw.res, GenoFile=genotypes, OutputFile=paste(out_base, "CFW.assoc", sep=""), control=marker_control)
GRAB.Marker(niger_seed_25.res, GenoFile=genotypes, OutputFile=paste(out_base, "Niger_Seed_25.assoc", sep=""), control=marker_control)
GRAB.Marker(niger_seed_30.res, GenoFile=genotypes, OutputFile=paste(out_base, "Niger_Seed_30.assoc", sep=""), control=marker_control)
GRAB.Marker(l_dopa_25.res, GenoFile=genotypes, OutputFile=paste(out_base, "L-DOPA_25.assoc", sep=""), control=marker_control)
GRAB.Marker(l_dopa_30.res, GenoFile=genotypes, OutputFile=paste(out_base, "L-DOPA-30.assoc", sep=""), control=marker_control)
GRAB.Marker(urease.res, GenoFile=genotypes, OutputFile=paste(out_base, "Urease.assoc", sep=""), control=marker_control)
GRAB.Marker(total_cell.res, GenoFile=genotypes, OutputFile=paste(out_base, "Total_Cell.assoc", sep=""), control=marker_control)
GRAB.Marker(cell_body.res, GenoFile=genotypes, OutputFile=paste(out_base, "Cell_Body.assoc", sep=""), control=marker_control)
GRAB.Marker(capsule.res, GenoFile=genotypes, OutputFile=paste(out_base, "Capsule.assoc", sep=""), control=marker_control)
GRAB.Marker(titan.res, GenoFile=genotypes, OutputFile=paste(out_base, "Titan.assoc", sep=""), control=marker_control)

