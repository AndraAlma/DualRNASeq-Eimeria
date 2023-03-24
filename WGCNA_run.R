install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library("WGCNA")
library("ComplexHeatmap")
library("edgeR")

### Get the file:
path_to_metadata_file <- "~/exjobb/full_metadata.csv"
path_to_chicken_reads <- "~/exjobb/chicken_counts.csv" 

# Read in count file 
chicken_reads <- read.delim(path_to_chicken_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)

metadata <- read.delim(path_to_metadata_file, sep = ";", check.names=FALSE, stringsAsFactors=FALSE)
chicken_metadata <- metadata 
rownames(chicken_metadata) <- chicken_metadata[,1]

# Read in count file 
chicken_reads <- read.delim(path_to_chicken_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)

# create groupings for timepoints
group_chicken <- factor(chicken_metadata[colnames(chicken_reads)[-c(1,2)],"Timepoint"])

#remove LOC112533599:
chicken_remove = which(chicken_reads[,1] == "LOC112533599") # Talk about why it's removed in the report, find other rRNA genes and remove them
chicken_reads = chicken_reads[-chicken_remove,]
#remove LOC112533601:
chicken_remove = which(chicken_reads[,1] == "LOC112533601") # Talk about why it's removed in the report, find other rRNA genes and remove them
chicken_reads = chicken_reads[-chicken_remove,]
# remove EYY68_mgr01
chicken_remove = which(chicken_reads[,1] == "EYY68_mgr01") # Talk about why it's removed in the report, find other rRNA genes and remove them
chicken_reads = chicken_reads[-chicken_remove,]
#Remove EYY68_mgr02
chicken_remove = which(chicken_reads[,1] == "EYY68_mgr02") # Talk about why it's removed in the report, find other rRNA genes and remove them
chicken_reads = chicken_reads[-chicken_remove,]

chicken_dgelist <- DGEList(counts = chicken_reads[,3:length(chicken_reads)], genes = chicken_reads[,c(1,2)], group = group_chicken)
rownames(chicken_dgelist$counts) <- rownames(chicken_dgelist$genes) <- chicken_dgelist$genes[,1]
keep_chicken <- filterByExpr(chicken_dgelist)
chicken_dgelist_filt <- chicken_dgelist[keep_chicken, ,] 
chicken_dgelist_filt$samples$lib.size <- colSums(chicken_dgelist_filt$counts)

# Normalize gene counts to account for a few genes dominating the counts of the samples
chicken_dgelist_filt_norm <- calcNormFactors(chicken_dgelist_filt)
chicken_cpm <- cpm.DGEList(chicken_dgelist_filt_norm)
chicken_cpm_df <- data.frame(chicken_cpm, check.names = FALSE)
m <- match(rownames(chicken_cpm_df), chicken_reads$gene_name)
chicken_cpm_df$entrez_gene_id <- chicken_reads$entrez_gene_id[m]
chicken_cpm_df <- chicken_cpm_df[,c(67,1:66)]

#option 1:tpm <- chicken_cpm_df
datExpr <- chicken_expr_cut

#simMatrix <- cor(t(tpm))
#datExpr=t(tpm[,2:67])

simMatrix <- cor(datExpr)
powers = c(c(1:10), seq(from = 11, to=20, by=1))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


plot_path <- "~/Exjobb/plots/WGCNA_soft_tresh.png"
png(plot_path, height = 600, width = 800, pointsize = 16)

par(mfrow=c(1,2))
#plot scale independance
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,col="red");
abline(h=0.80, col="red")
#Zoomed in
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"), ylim=range(0.6,0.95));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,col="red");
abline(h=0.80, col="red")

dev.off()

plot_path <- "~/Exjobb/plots/WGCNA_mean_con.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
#Connectivity
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"), ylim=range(0,500))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

dev.off()
### Compute adjacency matrices

beta <- 11
adjMatrix <- abs(simMatrix)^beta
unweightMatrix <- abs(simMatrix>=0.75)*1


### Plot Heatmaps of the PCC similarity matrix and the weighted and unweighted adjacency matrices derived from it.

s1 <- Heatmap(simMatrix[c(1:1000),c(1:1000)])
s2 <- Heatmap(adjMatrix[c(1:1000),c(1:1000)])
s3 <- Heatmap(unweightMatrix[c(1:1000),c(1:1000)])

hmlist <- s1 + s2 + s3
plot_path <- "~/Exjobb/plots/WGCNA_heatmaps_11.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
draw(hmlist)
dev.off()

TOM.w <- array(0, dim=c(0, ncol(adjMatrix), ncol(adjMatrix)))

### compute the TOM matrices
TOM.w <- TOMsimilarity(adjMatrix)

gc()
### clustering
wTree <- hclust(as.dist(1-TOM.w), method = "average")

### clustering using cutDynamicTree
w.clusters = cutreeDynamic(dendro = wTree, distM = 1-TOM.w, deepSplit = 2, cutHeight = 0.995, minClusterSize = 30, pamRespectsDendro = FALSE )
w.colors = labels2colors(w.clusters)

### summarize clusters
table(w.clusters)


### plot dendrograms and clusters
plotDendroAndColors(wTree, w.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

## Method 2:
options(stringsAsFactors = FALSE)

chicken_expr <- t(cpm(chicken_dgelist_filt_norm))

gsg = goodSamplesGenes(chicken_expr, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(chicken_expr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 50000, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)
keepSamples = (clust==1)
chicken_expr_cut = chicken_expr[keepSamples, ]
nGenes = ncol(chicken_expr_cut)
nSamples = nrow(chicken_expr_cut)

#powers = c(c(1:9), seq(from = 10, to=100, by=10))
# Call the network topology analysis function
#sft = pickSoftThreshold(chicken_expr_cut, powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2))
#cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
 #    ylim=range(0.25, 0.95), xlim=range(0,30), xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  #   main = paste("Scale independence"))
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    # labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
#abline(h=0.70,col="red")
#abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], sft$fitIndices[,5],
 #    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
 #    main = paste("Mean connectivity"))
#text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

