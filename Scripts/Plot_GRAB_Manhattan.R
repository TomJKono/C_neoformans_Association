# Make Manhattan plots of the GRAB results

setwd("/Users/tomkono/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/ST93A_Isolates/Manhattan_Plots/GRAB")

# Define the plotting characteristics
chrom_cols <- c("#000000", "#999999")
chrom_pch <- 19
chrom_cex <- 0.5

# Define the chromosome lengths - these will be approximations for conveneince
len_factor <- 1.02
chrom_lengths <- c(
    "1"=round(2205937*len_factor),
    "2"=round(1519281*len_factor),
    "3"=round(1334793*len_factor),
    "4"=round(1076465*len_factor),
    "5"=round(1781499*len_factor),
    "6"=round(1279254*len_factor),
    "7"=round(1341024*len_factor),
    "8"=round(1383765*len_factor),
    "9"=round(1115286*len_factor),
    "10"=round(989095*len_factor),
    "11"=round(1538971*len_factor),
    "12"=round(704606*len_factor),
    "13"=round(748576*len_factor),
    "14"=round(925085*len_factor)
    )
# And calculate a culumative position for plotting with offsets - these will
# be used for plotting the points.
chrom_lengths_cumulative <- c(0, cumsum(chrom_lengths)[1:length(chrom_lengths) - 1])
names(chrom_lengths_cumulative) <- names(chrom_lengths)
# Calculate cumulative positions again, but do not subtract off the length of
# the first chromosome. These will be used for drawing the labels on the
# X axis
chrom_lengths_cumulative_labs <- cumsum(chrom_lengths)
# Calculate a total length to define the plot bounds
total_length <- sum(chrom_lengths)
# Associate each chromosome name with a color:
chrom_col_vec <- sapply(
    seq_along(names(chrom_lengths)),
    function(x) {
        idx <- (x %% 2) + 1
        return(idx)}
        )
chrom_col_vec <- chrom_cols[chrom_col_vec]
names(chrom_col_vec) <- names(chrom_lengths)
# Calculate the midpoints of the chromosomes for labeling
chrom_midpoints <- sapply(
    names(chrom_lengths),
    function(x) {
        offset <- as.numeric(chrom_lengths_cumulative_labs[x])
        cmid <- as.numeric(chrom_lengths[x] / 2)
        return(offset - cmid)
    })

# Define a lil function to draw a plot from a GRAB P-values file.
manhattan_plot <- function(assoc_file) {
    gwas_res <- read.table(assoc_file, header=TRUE)
    pheno_name <- gsub(".assoc", "", basename(assoc_file))
    # Extract the position and the Wald P-value
    marker_info <- gwas_res$Info
    #   This column is formatted like CHROM:POS:REF:ALT
    #   We can use sapply() with the "[[" function to extract specific
    #   elements of the resulting list after calling strsplit(). Weird.
    chrom <- sapply(strsplit(marker_info, ":"), "[[", 1)
    pos <- as.numeric(sapply(strsplit(marker_info, ":"), "[[", 2))
    # pval <- -log10(gwas_res$Pvalue)
    pval <- gwas_res$Pvalue
    # Then calculate a FDR based on the P-value with the Benjamini-Hochberg
    # FDR prodedure
    fdr <- p.adjust(pval, method="BH")
    # Then we will take the log10 of 1/FDR
    l_i_fdr <- log10(1/fdr)
    # Calculate the cumulative position
    cumpos <- as.numeric(chrom_lengths_cumulative[chrom] + pos)
    # Choose the plotting color based on the chromosome name
    cols <- chrom_col_vec[chrom]
    # Draw the plot!
    #   Filename for the full dataset outputs
    # fname <- paste("Cn_38_Strains.", pheno_name, ".FDR.pdf", sep="")
    #   Filename for the ST93A subset outputs
    fname <- paste("ST93A.", pheno_name, ".FDR.pdf", sep="")
    pdf(file=fname, height=4, width=8)
    par(mar=c(4, 4, 2.5, 0.1))
    plot(0, xlab="", ylab="FDR", main=pheno_name, ylim=c(0.01, 2), xlim=c(0, total_length), type="n", axes=FALSE)
    points(
        l_i_fdr ~ cumpos,
        pch=chrom_pch,
        cex=chrom_cex,
        col=cols)
    axis(side=2, at=log10(1/c(1, 0.5, 0.25, 0.1, 0.05, 0.01)), labels=c(1, 0.5, 0.25, 0.1, 0.05, 0.01), las=2)
    axis(side=1, at=chrom_midpoints, labels=names(chrom_midpoints), cex.axis=0.8, las=2)
    # Add in dashed lines to separate the chromosomes
    abline(v=c(0, chrom_lengths_cumulative_labs), lwd=1, lty=3, col="darkgrey")
    abline(h=log10(1/0.05), lwd=1, lty=3, col="red")
    dev.off()
}

# And then apply this function over all the .assoc files!
#   This is for the full dataset
# assoc_dir <- "/Users/tomkono/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/Associations"
#   This is for the ST93A subset
assoc_dir <- "/Users/tomkono/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/ST93A_Isolates/Associations"
assoc_files <- list.files(assoc_dir, pattern="*.assoc$", full.names=TRUE)
dummy <- sapply(assoc_files, manhattan_plot)
