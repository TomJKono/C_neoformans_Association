# Make a pairs plot of the PCs as calculated by PLINK, a scree plot, and a
# labeled plot of the first two PCs (for talks etc). This version of the script
# will produce plots using all variants, rather than just the 692 that are in
# gene coding regions.

library(ggplot2)
library(ggrepel)

# For HST laptop
setwd("/Users/konox006/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/PCA")
# For non-HST laptop
# setwd("/Users/tomkono/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/ST93A_Isolates/PCA")

#   These are for the full dataset
eigvec <- read.table("../../../Data/38_Strain_PLINK/MAF_Filtered/All_Variants_PCA_woOutgroup.eigenvec", header=FALSE)
eigval <- read.table("../../../Data/38_Strain_PLINK/MAF_Filtered/All_Variants_PCA_woOutgroup.eigenval", header=FALSE)
#   These are for the ST93A subset
# eigvec <- read.table("../../../Data/38_Strain_PLINK/ST93A_Isolates/ST93A_Isolates.PCA.eigenvec", header=FALSE)
# eigval <- read.table("../../../Data/38_Strain_PLINK/ST93A_Isolates/ST93A_Isolates.PCA.eigenval", header=FALSE)

# The first two columns are family ID and individual ID; number of PCs is the
# total column number minus 2
npcs <- ncol(eigvec) - 2
colnames(eigvec) <- c("FID", "IID", paste("PC", seq(1, npcs), sep=""))

# pdf(file="Cn_38_Strains.AllVariants.ST93A.PCA.pdf", height=12, width=12)
pdf(file="Cn_38_Strains.AllVariants.PCA.woOutgroup.pdf", height=12, width=12)
pairs(
    eigvec[, -c(1, 2)],
    upper.panel=NULL,
    pch=19,
    cex=0.75)
dev.off()

# Calcluate the proportion of variance explained by each principal component.
# Based on this, we will use the first two PCs as covariates.
# pdf(file="Cn_38_Strains.AllVariants.ST93A.Scree.woOutgroup.pdf", height=6, width=6)
pdf(file="Cn_38_Strains.AllVariants.Scree.woOutgroup.pdf", height=6, width=6)
pvar <- eigval$V1 / sum(eigval$V1)
labs <- paste("PC", seq(1, npcs), sep="")
at <- barplot(
    pvar,
    xlab="Principal Component",
    ylab="Proportion of Variance Explained",
    main="38 Strain Scree Plot",
    axes=FALSE)
axis(side=2)
axis(side=1, at=at, labels=labs, las=2, cex.axis=0.75)
dev.off()

# Then, we will draw a labeled plot of PC1 vs PC2
# pdf(file="Cn_38_Strains.AllVariants.ST93A.PC1-vs-PC2.woOutgroup.pdf", height=6, width=6)
pdf(file="Cn_38_Strains.AllVariants.PC1-vs-PC2.woOutgroup.pdf", height=6, width=6)
plt <- ggplot(eigvec, aes(x=PC1, y=PC2, label=IID)) +
    geom_point(color="red") +
    xlim(-0.3, 0.4) +
    ylim(-0.25, 1) +
    xlab(paste("PC1 (", round(pvar[1]*100, 2), "% Variance Explained)", sep="")) +
    ylab(paste("PC2 (", round(pvar[2]*100, 2), "% Variance Explained)", sep="")) +
    geom_text_repel(
        box.padding=0.3,
        max.overlaps=Inf,
        segment.size=0.1,
        force=1,
        force_pull=0.1,
        size=2.5) +
    theme_bw()
plt
dev.off()
