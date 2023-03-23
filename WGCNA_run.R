install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library("WGCNA")
library("ComplexHeatmap")
tpm <- chicken_cpm_df

simMatrix <- cor(t(tpm))
datExpr=t(tpm[,2:67])
powers = c(c(2:10), seq(from = 11, to=20, by=1))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

par(mfrow=c(1,2))

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,col="red");
abline(h=0.80, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

### Load libraries

library(ComplexHeatMap)
### Compute adjacency matrices

beta <- 11
adjMatrix <- abs(simMatrix)^beta
unweightMatrix <- abs(simMatrix>=0.75)*1


### Plot Heatmaps of the PCC similarity matrix and the weighted and unweighted adjacency matrices derived from it.

s1 <- Heatmap(simMatrix[c(1:1000),c(1:1000)])
s2 <- Heatmap(adjMatrix[c(1:1000),c(1:1000)])
s3 <- Heatmap(unweightMatrix[c(1:1000),c(1:1000)])

hmlist <- s1 + s2 + s3
plot_path <- "~/Exjobb/plots/WGCNA_heatmaps.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
draw(hmlist)
dev.off()

TOM.w <- array(0, dim=c(0, ncol(adjMatrix), ncol(adjMatrix)))

### compute the TOM matrices
TOM.w <- TOMsimilarity(adjMatrix)

### clustering
wTree <- hclust(as.dist(1-TOM.w), method = "average")

### clustering using cutDynamicTree
w.clusters = cutreeDynamic(dendro = wTree, distM = 1-TOM.w, deepSplit = 2, cutHeight = 0.995, minClusterSize = 30, pamRespectsDendro = FALSE )


w.colors = labels2colors(w.clusters)


### summarize clusters


table(w.clusters)


w.clusters


0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 


7 4766 4358 1955 1952 1144 744 619 603 559 512 498 159 115 92 88 


16 17 


85 61 


### plot dendrograms and clusters


plotDendroAndColors(wTree, w.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

