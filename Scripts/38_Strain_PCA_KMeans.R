# Run a K-means clustering on the strains, using the 2 PCs as the
# observations. We are doing this just to show that the inference of 3
# clusters of data is statistically sound. We will use the 692 variants because
# the clustering inference was made off of those.

library(cluster)
library(factoextra)
library(ggrepel)

set.seed(54321)

# For HST laptop
# setwd("/Users/konox006/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/PCA")
# For non-HST laptop
setwd("/Users/tomkono/Dropbox/GitHub/RIS/Nielsen_Cn_GWAS/Results/38_Strain_GWAS/PCA")

# Read in the PC coordinates
eigvec <- read.table("../../../Data/38_Strain_PLINK/MAF_Filtered/Cn_38_Strains.MAF.PCA.eigenvec", header=FALSE)
eigval <- read.table("../../../Data/38_Strain_PLINK/MAF_Filtered/Cn_38_Strains.MAF.PCA.eigenval", header=FALSE)

# The first two columns are family ID and individual ID; number of PCs is the
# total column number minus 2
npcs <- ncol(eigvec) - 2
colnames(eigvec) <- c("FID", "IID", paste("PC", seq(1, npcs), sep=""))

# Then, we want to make a dataframe where the rows are strains and the columns
# are the PC values. We are close:
clust_in <- as.data.frame(eigvec[, -c(1:2, 5:12)])
rownames(clust_in) <- eigvec$IID
print(clust_in)

# Let's visualize Euclidean distances between strains
euc_dist <- get_dist(clust_in, method="euclidean")
pdf(file="PC_Euclidean_Distances.pdf", height=8, width=8)
fviz_dist(euc_dist)
dev.off()

# Now, let's define a range of k-means across which to cluster and calculate
# the "gap statistic"
gap_stats <- clusGap(clust_in, FUN=kmeans, nstart=50, K.max=10, B=1000)
pdf(file="Gap_Statistic.pdf", height=6, width=6)
fviz_gap_stat(gap_stats)
dev.off()

# Let's calculate the silhouette score, too
avg_silhouette <- function(kval, clust_dat) {
    k_clust <- kmeans(clust_dat, centers=kval, nstart=50)
    s_dist <- silhouette(k_clust$cluster, dist(clust_dat))
    return(mean(s_dist[,3]))
}
k_vals <- 2:10
s_dists <- sapply(k_vals, avg_silhouette, clust_in)
pdf(file="Silhouette_Scores.pdf", height=6, width=6)
plot(
    s_dists ~ k_vals,
    xlim=c(1, 10),
    ylim=c(-1, 1),
    xlab="K-value",
    ylab="Average Silhouette Score",
    main="K-means Clustering of PC Distance",
    pch=19)
abline(h=0, lwd=1, lty=3, col="black")
abline(v=3, lwd=1, lty=3, col="blue")
dev.off()

# And finally, let's plot the clusters with labels
fin_clust <- kmeans(clust_in, centers=3, nstart=50, iter.max=1000)
clustered <- clust_in
clustered$IID <- rownames(clust_in)
clustered$Cluster <- factor(fin_clust$cluster)
pvar <- eigval$V1 / sum(eigval$V1)
pdf(file="KMeans_PCA_Labeled.pdf", height=6, width=6)
plt <- ggplot(clustered, aes(x=PC1, y=PC2, label=IID, color=Cluster)) +
    geom_point() +
    xlim(-0.4, 0.7) +
    ylim(-0.5, 0.5) +
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
