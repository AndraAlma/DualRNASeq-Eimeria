install.packages("BiocManager")
BiocManager::install("ggbreak")

#Load in correct packages:

library("DESeq2")
library("ggbiplot")
library("ggplot2")
library("edgeR")
library("statmod")
library("EnhancedVolcano")
library("Rmisc")
library("org.Gg.eg.db")
library("ReportingTools")
library("ClassDiscovery")
library("RColorBrewer")
library("gplots")
library("devtools")
library("KEGGREST")
library("WGCNA")
library("ComplexHeatmap")
library("magrittr")
library("ggforce")
library("readr")
library("ggbreak")
install_github("vqv/ggbiplot")


# File paths to count & metadata:
path_to_metadata_file <- "full_metadata.csv"
path_to_chicken_reads <- "chicken_counts.csv" 


# Read in files:
metadata <- read.delim(path_to_metadata_file, sep = ";", check.names=FALSE, stringsAsFactors=FALSE)
chicken_metadata <- metadata 
rownames(chicken_metadata) <- chicken_metadata[,1]
chicken_reads <- read.delim(path_to_chicken_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)


# Create groupings for timepoints:
group_chicken <- factor(chicken_metadata[colnames(chicken_reads)[-c(1,2)],"Timepoint"])


# Normalize counts to counts per million:
chicken_cpm_raw <- cpm(chicken_reads[,-c(1,2)])
rownames(chicken_cpm_raw) <- chicken_reads[,1]


# Check PCA of normalized counts to detect outliers:
pca_chicken <- prcomp(chicken_cpm_raw)
ggbiplot(pca_chicken, var.axes = FALSE) + coord_cartesian(xlim = c(0, 160)) +
  geom_text(aes(label = chicken_reads[,1]), hjust = 0.5, vjust = -0.5, size = 3) +
  labs(title = "PCA of normalized chicken read counts per gene") +
  theme(plot.title = element_text(hjust = 0.5))


# Remove outliers:
chicken_remove = which(chicken_reads[,1] == "LOC112533599")
chicken_reads = chicken_reads[-chicken_remove,]

chicken_remove = which(chicken_reads[,1] == "LOC112533601") 
chicken_reads = chicken_reads[-chicken_remove,]

chicken_remove = which(chicken_reads[,1] == "EYY68_mgr01")
chicken_reads = chicken_reads[-chicken_remove,]

chicken_remove = which(chicken_reads[,1] == "EYY68_mgr02")
chicken_reads = chicken_reads[-chicken_remove,]

chicken_cpm_raw <- cpm(chicken_reads[,-c(1,2)])
rownames(chicken_cpm_raw) <- chicken_reads[,1]


# Check PCA of normalized counts again:
pca_chicken <- prcomp(chicken_cpm_raw)
ggbiplot(pca_chicken, var.axes = FALSE) + coord_cartesian(xlim = c(0, 80)) +
  geom_text(aes(label = chicken_reads[,1]), hjust = 0.5, vjust = -0.5, size = 3) +
  labs(title = "PCA of normalized chicken read counts per gene") +
  theme(plot.title = element_text(hjust = 0.5))


# Filter out low-expressed genes:
chicken_dgelist <- DGEList(counts = chicken_reads[,3:length(chicken_reads)], genes = chicken_reads[,c(1,2)], group = group_chicken)
rownames(chicken_dgelist$counts) <- rownames(chicken_dgelist$genes) <- chicken_dgelist$genes[,1]
keep_chicken <- filterByExpr(chicken_dgelist)
chicken_dgelist_filt <- chicken_dgelist[keep_chicken, ,] 
chicken_dgelist_filt$samples$lib.size <- colSums(chicken_dgelist_filt$counts)


# Normalize gene counts to account for a few genes dominating the counts of the samples:
chicken_dgelist_filt_norm <- calcNormFactors(chicken_dgelist_filt)
chicken_labels <- paste(rownames(chicken_dgelist_filt_norm$samples), chicken_dgelist_filt_norm$samples$group, sep = "_")


# Plot multidimentional scaling:
mds_path <- "plots/mds_plots.png"
png(mds_path, height = 1200, width = 800, pointsize = 16)
par(mfrow=c(2,1))
p_chicken <- plotMDS(chicken_dgelist_filt_norm, 
                     labels = chicken_dgelist_filt_norm$samples$group)
dev.off()

chicken_mds_df <- data.frame(x = p_chicken$x, 
                             y = p_chicken$y, 
                             labels = as.character(chicken_dgelist_filt_norm$samples$group), 
                             stringsAsFactors = FALSE)

p_chicken_gg <- ggplot(chicken_mds_df, aes(x = x, y = y, color=labels)) +
  geom_point(size = 5) +
  scale_fill_brewer(palette="Dark2")+
  coord_cartesian(ylim = c(-1.8, 1.3)) +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.position = "right")

mds_path <- "plots/mds_plots.pdf"
pdf(mds_path, height = 12, width = 12)
p_chicken_gg
dev.off()


# Normalize to counts per million, add entrez gene IDs and print normalized counts to table:
chicken_cpm <- cpm.DGEList(chicken_dgelist_filt_norm)
chicken_cpm_df <- data.frame(chicken_cpm, check.names = FALSE)
m <- match(rownames(chicken_cpm_df), chicken_reads$gene_name)
chicken_cpm_df$entrez_gene_id <- chicken_reads$entrez_gene_id[m]
chicken_cpm_df <- chicken_cpm_df[,c(67,1:66)]

write.csv(chicken_cpm_df, file = "tables/chicken_counts_norm.csv")


# Plot PCA:
pca_cpm_chicken <- prcomp(t(chicken_cpm))
pca_path <- "plots/chicken_pca.png"
png(pca_path, height = 500, width = 520)
ggbiplot(pca_cpm_chicken, var.axes = FALSE, choices = 1:2, alpha = 1) +
  geom_text(aes(label = chicken_dgelist_filt_norm$samples$group), hjust = 0.5, vjust = -0.5, size = 4.5) +
  ggtitle("PCA of normalized chicken read counts per sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()


# Create design matrix:
chicken_dgelist_filt_norm$samples$group <- relevel(chicken_dgelist_filt_norm$samples$group, ref="0")
design_chicken <- model.matrix(~0+group, data = chicken_dgelist_filt_norm$samples)
rownames(design_chicken) <- colnames(chicken_dgelist_filt_norm)

chicken_disp <- estimateDisp(chicken_dgelist_filt_norm, design_chicken, robust = TRUE)
plotBCV(chicken_disp)


# Create contrasts for DE:
fit_chicken <- glmQLFit(chicken_disp, design_chicken, prior.count = 0.125)
contrasts_chicken <- makeContrasts(Uvs1=group1-group0, Uvs2=group2-group0, Uvs3=group3-group0, Uvs4=group4-group0, Uvs10=group10-group0,
                                   levels = design_chicken)

qlf_chicken_Uvs1 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs1"])
qlf_chicken_Uvs2 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs2"])
qlf_chicken_Uvs3 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs3"])
qlf_chicken_Uvs4 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs4"])
qlf_chicken_Uvs10 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs10"])

chicken_qlf_list <- list(qlf_chicken_Uvs1,qlf_chicken_Uvs2,qlf_chicken_Uvs3,qlf_chicken_Uvs4,
                         qlf_chicken_Uvs10)
topgenes_chicken_list <- lapply(chicken_qlf_list,
                                function(x) topTags(x, n = dim(chicken_qlf_list[[1]]$table)[1]))
 

# Define FDR & logFC threshold for significance:
fdr_threshold <- 0.05
logfc_threshold <- 1


# Create volcano plots for all comparisons:
timepoints <- c("1_day","2_days","3_days","4_days","10_days")
plot_list <- list()
i = 1
while (i <= length(topgenes_chicken_list)) {
  volcano_path <- paste("plots/chicken_volcano_lab_", timepoints[i], ".pdf", sep = "")
  pdf(volcano_path, height = 10, width = 16)
  p <- EnhancedVolcano(topgenes_chicken_list[[i]]$table,
                       lab = rownames(topgenes_chicken_list[[i]]$table),
                       title = timepoints[i],
                       subtitle = "",
                       x = 'logFC',
                       y = 'FDR',
                       xlim = c(-4,4),
                       ylim = c(0,7),
                       ylab = bquote(~-Log[10]~italic(FDR)),
                       legendLabels = c('NS', expression(Log[2]~FC),
                                        "FDR", expression(FDR~and~log[2]~FC)),
                       caption = "",
                       pointSize = 4,
                       labSize = 5,
                       pCutoff = fdr_threshold,
                       FCcutoff = logfc_threshold) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 16))
  plot_list[[i]] <- p
  print(p)
  dev.off()
  i = i + 1
}
volcano_path <- "plots/all_chicken_plots_lab.png"
png(volcano_path, height = 1200, width = 1200)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

volcano_path <- "plots/all_chicken_plots_lab.pdf"
pdf(volcano_path, height = 20, width = 20)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()


# Volcano plots without labels:
i = 1
while (i <= length(topgenes_chicken_list)) {
  volcano_path <- paste("plots/chicken_volcano_", timepoints[i], ".pdf", sep = "")
  pdf(volcano_path, height = 10, width = 16)
  p <- EnhancedVolcano(topgenes_chicken_list[[i]]$table,
                       lab = rownames(topgenes_chicken_list[[i]]$table),
                       title = timepoints[i],
                       subtitle = "",
                       x = 'logFC',
                       y = 'FDR',
                       xlim = c(-4,4),
                       ylim = c(0,7),
                       ylab = bquote(~-Log[10]~italic(FDR)),
                       legendLabels = c('NS', expression(Log[2]~FC),
                                        "FDR", expression(FDR~and~log[2]~FC)),
                       caption = "",
                       pointSize = 4,
                       labSize = 0,
                       pCutoff = fdr_threshold,
                       FCcutoff = logfc_threshold) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 16))
  plot_list[[i]] <- p
  print(p)
  dev.off()
  i = i + 1
}
volcano_path <- "plots/all_chicken_plots.png"
png(volcano_path, height = 1200, width = 1200)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

volcano_path <- "plots/all_chicken_plots.pdf"
pdf(volcano_path, height = 20, width = 20)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()


# Compute the number of unique genes differentially expressed at each time point, 
# which genes they are and which genes are differentially expressed at all time points:
get_unique_and_shared_de_genes <- function(de_gene_list, id_col_num = 1) {
  unique_de_gene_list <- de_gene_list
  shared_genes <- c()
  num_shared_genes <- 0
  
  i <- 0
  while (i < length(de_gene_list)) {
    i <- i + 1
    k <- 0
    while (k < dim(unique_de_gene_list[[i]])[1]) {
      j <- 0
      num_samples <- 0
      while (j < length(de_gene_list)) {
        j <- j + 1
        if (j == i) {
          next()
        } else {
          name_in_k <- unique_de_gene_list[[i]][k+1,id_col_num] %in% de_gene_list[[j]][,id_col_num]
        }
        if (name_in_k) {
          num_samples <- num_samples + 1
        }
      }
      if (num_samples == 5 && !(unique_de_gene_list[[i]][k+1,id_col_num] %in% shared_genes)) {
        num_shared_genes <- num_shared_genes + 1
        shared_genes[num_shared_genes] <- unique_de_gene_list[[i]][k+1,id_col_num]
      }
      if (num_samples > 0) {
        unique_de_gene_list[[i]] <- unique_de_gene_list[[i]][-(k+1),]
      } else {
        k <- k + 1
      }
    }
  }
  return(list(unique_de_gene_list, shared_genes))
}


de_genes_chicken <- lapply(topgenes_chicken_list, function(x) x$table[x$table$FDR < fdr_threshold,])
de_genes_chicken_pos <- lapply(de_genes_chicken, function(x) x[x$logFC > logfc_threshold,])
de_genes_chicken_neg <- lapply(de_genes_chicken, function(x) x[x$logFC < -logfc_threshold,])
de_genes_chicken_all <- lapply(de_genes_chicken, function(x) x[abs(x$logFC) > logfc_threshold,])

unique_chicken_de_genes_pos <- get_unique_and_shared_de_genes(de_genes_chicken_pos, 1)
unique_chicken_de_genes_neg <- get_unique_and_shared_de_genes(de_genes_chicken_neg, 1)
unique_chicken_de_genes_all <- get_unique_and_shared_de_genes(de_genes_chicken_all, 1)


# Get a list of all genes that are significantly DE at any timepoint and a list 
# of those that are DE and have a logFC > 1:
egGENENAME <- toTable(org.Gg.egGENENAME)
i <- 0
while (i < length(de_genes_chicken)) {
  i <- i + 1
  df_temp <- de_genes_chicken[[i]][,c("gene_name", "entrez_gene_id", "FDR")]
  if (i == 1) {
    all_de_genes_chicken <- df_temp
  } else {
    all_de_genes_chicken <- rbind(all_de_genes_chicken, df_temp)
  }
}
all_de_genes_chicken <- all_de_genes_chicken[order(all_de_genes_chicken$FDR),]
all_de_genes_chicken <- all_de_genes_chicken[!duplicated(all_de_genes_chicken[,c("gene_name","entrez_gene_id")]),]
all_de_genes_chicken$description <- egGENENAME[match(all_de_genes_chicken$entrez_gene_id, egGENENAME$gene_id),]$gene_name

write.table(all_de_genes_chicken, "tables/all_de_genes_chicken.csv", sep = ",", row.names = FALSE)

i <- 0
while (i < length(de_genes_chicken_all)) {
  i <- i + 1
  df_temp <- de_genes_chicken_all[[i]][,c("gene_name", "entrez_gene_id", "FDR")]
  if (i == 1) {
    logfc_filt_de_genes_chicken <- df_temp
  } else {
    logfc_filt_de_genes_chicken <- rbind(logfc_filt_de_genes_chicken, df_temp)
  }
}
logfc_filt_de_genes_chicken <- logfc_filt_de_genes_chicken[order(logfc_filt_de_genes_chicken$FDR),]
logfc_filt_de_genes_chicken <- logfc_filt_de_genes_chicken[!duplicated(logfc_filt_de_genes_chicken[,c("gene_name","entrez_gene_id")]),]
rownames(logfc_filt_de_genes_chicken) <- logfc_filt_de_genes_chicken$gene_name
logfc_filt_de_genes_chicken$description <- egGENENAME[match(logfc_filt_de_genes_chicken$entrez_gene_id, egGENENAME$gene_id),]$gene_name
missing <- which(is.na(logfc_filt_de_genes_chicken$description))
i <- 0
while (i <= length(missing)) {
  logfc_filt_de_genes_chicken$description[missing[i]]<-"NA"
  i <- i+1
}
write.table(logfc_filt_de_genes_chicken, "tables/logfc_filt_de_genes_chicken.csv", sep = ",", row.names = FALSE)


# Get the top 30 most significantly DE genes from each timepoint and export to a table:
top_de_genes_path <- "tables"
top_de_gene_file_name <- "top_de_genes"
logfc_filtered_topgenes_chicken_list <- lapply(topgenes_chicken_list, function(x) x[abs(x$table$logFC) > logfc_threshold, ])

i <- 0
while (i < length(timepoints)) {
  i <- i + 1
  top_de_genes_file_path_chicken <- paste(top_de_genes_path, "/", top_de_gene_file_name, "_chicken_", timepoints[i], ".csv", sep = "")
  write.table(logfc_filtered_topgenes_chicken_list[[i]][1:30,], top_de_genes_file_path_chicken, sep = ",", row.names = FALSE)
}


# Filter for FDR as well, and print top 40 to a list:
top_de_genes_path <- "tables"
top_de_gene_file_name <- "top_de_genes_FDR"
i <- 0
while (i < length(timepoints)) {
  i <- i + 1
  top_de_genes_file_path_chicken <- paste(top_de_genes_path, "/", top_de_gene_file_name, "_chicken_", timepoints[i], ".csv", sep = "")
  write.table(de_genes_chicken_all[[i]][1:40,], top_de_genes_file_path_chicken, sep = ",", row.names = FALSE)
}


# Set up functions to create heatmaps:
create_gene_heatmap_df <- function(de_gene_list, heatmap_value, gene_names, timepoints) {
  # de_gene_list is a list of DEGlist objects after performing a DE analysis
  # heatmap_value is a string with a value of either "FDR" or "logFC" depending on what the heatmap should show
  # gene_names is a vector of gene names in the analysis
  # timepoints is a vector of time points at which the samples were taken
  mapped_val_df <- as.data.frame(matrix(0,
                                        ncol = length(de_gene_list),
                                        nrow = dim(de_gene_list[[1]])[1]))
  
  colnames(mapped_val_df) <- timepoints
  
  rownames(mapped_val_df) <- gene_names
  if (heatmap_value == "FDR") {
    i <- 0
    while (i < length(de_gene_list)) {
      i <- i + 1
      m <- match(gene_names, rownames(de_gene_list[[i]]$table))
      mapped_val_df[,i] <- -log10(de_gene_list[[i]]$table[m,]$FDR)
    }
  } else if(heatmap_value == "logFC") {
    i <- 0
    while (i < length(de_gene_list)) {
      i <- i + 1
      m <- match(gene_names, rownames(de_gene_list[[i]]$table))
      mapped_val_df[,i] <- de_gene_list[[i]]$table[m,]$logFC
    }
  }
  return(mapped_val_df)
}

filter_gene_heatmap_df <- function(gene_heatmap_df, de_gene_list) {
  m <- match(de_gene_list, rownames(gene_heatmap_df))
  m <- m[complete.cases(m)]
  return(gene_heatmap_df[m,])
}


# Create heatmaps for all genes based on logFC and FDR:
chicken_fdr_df <- create_gene_heatmap_df(topgenes_chicken_list, "FDR", rownames(topgenes_chicken_list[[1]]$table), timepoints)
chicken_logfc_df <- create_gene_heatmap_df(topgenes_chicken_list, "logFC", rownames(topgenes_chicken_list[[1]]$table), timepoints)
chicken_fdr_df_filt <- filter_gene_heatmap_df(chicken_fdr_df, rownames(logfc_filt_de_genes_chicken))
chicken_logfc_df_filt <- filter_gene_heatmap_df(chicken_logfc_df, rownames(logfc_filt_de_genes_chicken))

heatmap_path <- "plots/chicken_genes_pval.pdf"
pdf(heatmap_path, width = 30, height = 240, pointsize = 60)
aspectHeatmap(as.matrix(chicken_fdr_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene p-values", 
              col = colorRampPalette(brewer.pal(9, "YlOrRd"))(500),
              hExp = 1, wExp = 1)
dev.off()

heatmap_path <- "plots/chicken_genes_pval.png"
png(heatmap_path, height = 1200, width = 1200)
aspectHeatmap(as.matrix(chicken_fdr_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene p-values", 
              col = colorRampPalette(brewer.pal(9, "YlOrRd"))(500),
              hExp = 1, wExp = 1)
dev.off()

heatmap_path <- "plots/chicken_genes_logfc.pdf"
pdf(heatmap_path, width = 240, height = 180, pointsize = 100)
aspectHeatmap(as.matrix(chicken_logfc_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene log fold change", 
              col = colorRampPalette(brewer.pal(9, "RdYlGn"))(500),
              hExp = 1, wExp = 1)
dev.off()

heatmap_path <- "plots/chicken_genes_logfc.png"
png(heatmap_path, height = 1200, width = 1200)
aspectHeatmap(as.matrix(chicken_logfc_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene log fold change", 
              col = colorRampPalette(brewer.pal(9, "RdYlGn"))(500),
              hExp = 1, wExp = 1)
dev.off()


# GO category and KEGG pathway analysis of the data.  
# For each list, the data is in order of increasing time
go_chicken_list <- lapply(chicken_qlf_list, 
                          function(x) goana(x, geneid = chicken_dgelist_filt_norm$genes$entrez_gene_id, FDR = fdr_threshold, species = "Gg"))
kegg_chicken_list <- lapply(chicken_qlf_list, 
                            function(x) kegga(x, geneid = chicken_dgelist_filt_norm$genes$entrez_gene_id, FDR = fdr_threshold, species.KEGG = "gga"))

topgo_chicken_up_list <- lapply(go_chicken_list, 
                                function(x) topGO(x, ont="BP", sort="Up", n=250))
topgo_chicken_down_list <- lapply(go_chicken_list, 
                                  function(x) topGO(x, ont="BP", sort="Down", n=250))
topgo_chicken_all_list <- lapply(go_chicken_list, 
                                 function(x) topGO(x, ont="BP", n=250))
topkegg_chicken_up_list <- lapply(kegg_chicken_list,
                                  function(x) topKEGG(x, sort="Up", n=250))
topkegg_chicken_down_list <- lapply(kegg_chicken_list,
                                    function(x) topKEGG(x, sort="Down", n=250))
topkegg_chicken_all_list <- lapply(kegg_chicken_list, 
                                   function(x) topKEGG(x, n=250))


# Print a list of the top 50 lowest p-value categories for each time point
i <- 0
while (i < length(topgo_chicken_all_list)) {
  i <- i + 1
  go_path <- paste("tables/top_go_", timepoints[i], ".csv", sep = "")
  kegg_path <- paste("tables/top_kegg_", timepoints[i], ".csv", sep = "")
  write.csv(topgo_chicken_all_list[[i]][1:50,], go_path)
  write.csv(topkegg_chicken_all_list[[i]][1:50,], kegg_path)
}


# Set up functions to create GO and KEGG category heatmaps:
create_cat_pval_df <- function(cat_list, cat_name, timepoints) {
  # cat_name_col is a vector containing names of categories
  cat_pval_up_df <- as.data.frame(matrix(0,
                                         ncol = length(cat_list),
                                         nrow = dim(cat_list[[1]])[1]))
  cat_pval_down_df <- as.data.frame(matrix(0,
                                           ncol = length(cat_list),
                                           nrow = dim(cat_list[[1]])[1]))
  cat_pval_mix_df <- as.data.frame(matrix(0,
                                          ncol = length(cat_list),
                                          nrow = dim(cat_list[[1]])[1]))
  colnames(cat_pval_up_df) <- timepoints
  colnames(cat_pval_down_df) <- timepoints
  colnames(cat_pval_mix_df) <- timepoints
  
  rownames(cat_pval_up_df) <- cat_name
  rownames(cat_pval_down_df) <- cat_name
  rownames(cat_pval_mix_df) <- cat_name
  
  i <- 0
  while (i < length(cat_list)) {
    i <- i + 1
    cat_pval_up_df[,i] <- -log10(cat_list[[i]]$P.Up)
    cat_pval_down_df[,i] <- -log10(cat_list[[i]]$P.Down)
  }
  
  i <- 0
  while (i < dim(cat_pval_up_df)[1]) {
    i <- i + 1
    j <- 0
    while (j < dim(cat_pval_up_df)[2]) {
      j <- j + 1
      if (cat_pval_up_df[i,j] >= cat_pval_down_df[i,j]) {
        cat_pval_mix_df[i,j] <- cat_pval_up_df[i,j]
      } else {
        cat_pval_mix_df[i,j] <- -cat_pval_down_df[i,j]
      }
    }
  }
  return(list(cat_pval_up_df, cat_pval_down_df, cat_pval_mix_df))
}

remove_high_pval_rows <- function(pval_df, pval_thresh) {
  pval_thresh <- -log10(pval_thresh)
  i <- 1
  while (i <= dim(pval_df)[1]) {
    num_high_pval <- 0
    j <- 0
    while (j < dim(pval_df)[2]) {
      j <- j + 1
      if (pval_df[i,j] < pval_thresh) {
        num_high_pval <- num_high_pval + 1
      }
    }
    if (num_high_pval == j) {
      pval_df <- pval_df[-i,]
    } else {
      i <- i + 1
    }
  }
  return(pval_df)
}


# Create heatmaps for GO and KEGG categories:
go_chicken_pval_list <- create_cat_pval_df(go_chicken_list, go_chicken_list[[1]]$Term, timepoints)
kegg_chicken_pval_list <- create_cat_pval_df(kegg_chicken_list, kegg_chicken_list[[1]]$Pathway, timepoints)

go_chicken_pval_list_filtered <- lapply(go_chicken_pval_list, function(x) remove_high_pval_rows(x, fdr_threshold/10))
kegg_chicken_pval_list_filtered <- lapply(kegg_chicken_pval_list, function(x) remove_high_pval_rows(x, fdr_threshold))

heatmap_path <- "plots/GO_cat_upregulated.pdf"
pdf(heatmap_path, width = 60, height = 180, pointsize = 60)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated GO categories", 
              hExp = 5, wExp = 0.8, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "plots/GO_cat_downregulated.pdf"
pdf(heatmap_path, width = 140, height = 180, pointsize = 100)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[2]]), Colv = NA,
              xlab = "Timepoints", main = "Downregulated GO categories", 
              hExp = 2, wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-3,length = 15), seq(-2.95, 2.95, length = 119), seq(3,10,length = 15))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 148)

heatmap_path <- "plots/GO_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 60)
heatmap.2(as.matrix(go_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
          xlab = "Timepoints", main = "GO categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45)
dev.off()

heatmap_path <- "plots/KEGG_cat_upregulated.pdf"
pdf(heatmap_path, width = 130, height = 120, pointsize = 120)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated KEGG categories", 
              wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "plots/KEGG_cat_downregulated.pdf"
pdf(heatmap_path, width = 100, height = 100, pointsize = 100)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[2]]), Colv = NA,
              xlab = "Timepoints", main = "Downregulated KEGG categories", 
              hExp = 2, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-5,length = 11), seq(-4.95, 4.95, length = 199), seq(5,10,length = 11))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 220)

heatmap_path <- "plots/KEGG_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 100)
heatmap.2(as.matrix(kegg_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
          xlab = "Timepoints", main = "KEGG categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45)
dev.off()



# WGCNA gene network analysis

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



# Plot sample tree:
sampleTree = hclust(dist(chicken_expr), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 50000, col = "red")


# Determine cluster under the treshold line
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)
keepSamples = (clust==1)
chicken_expr_cut = chicken_expr[keepSamples, ]
nGenes = ncol(chicken_expr_cut)
nSamples = nrow(chicken_expr_cut)

datExpr <- chicken_expr_cut
simMatrix <- cor(datExpr)
names(datExpr) = colnames(chicken_expr_cut);
rownames(datExpr) = rownames(chicken_expr_cut)
rm(chicken_expr_cut)
rm(chicken_cpm)

# Find soft treshold:
powers = c(c(1:10), seq(from = 11, to=20, by=1))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


# Plot soft treshold scale independence:
plot_path <- "plots/WGCNA_soft_tresh.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,col="red");
abline(h=0.80, col="red")
#Zoomed in
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"), ylim=range(0.6,0.95));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,col="red");
abline(h=0.80, col="red")
dev.off()


# Plot mean connectivity:
plot_path <- "plots/WGCNA_mean_con.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")
# Zoomed in
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"), ylim=range(0,500))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")
dev.off()


# Compute adjacency matrices:
beta <- 11
adjMatrix <- abs(simMatrix)^beta
unweightMatrix <- abs(simMatrix>=0.75)*1


# Plot Heatmaps of the PCC similarity matrix and the weighted and unweighted adjacency matrices derived from it.
s1 <- Heatmap(simMatrix[c(1:1000),c(1:1000)])
s2 <- Heatmap(adjMatrix[c(1:1000),c(1:1000)])
s3 <- Heatmap(unweightMatrix[c(1:1000),c(1:1000)])

hmlist <- s1 + s2 + s3
plot_path <- "plots/WGCNA_heatmaps_11.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
draw(hmlist)
dev.off()

rm(unweightMatrix)
rm(simMatrix)
rm(adjMatrix)

cor <- WGCNA::cor
# Run blockwise modules:
net11 <- blockwiseModules(datExpr, power = 11,
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "fusionTOM9", 
                          verbose = 3,
                          maxBlockSize = 17000)



table(net11$colors)

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net11$colors)
# Plot the dendrogram and the module colors underneath
plot_path <- "plots/WGCNA_modules_block_net.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
plotDendroAndColors(net11$dendrograms[[1]], mergedColors[net11$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net11$colors
moduleColors = labels2colors(net11$colors)

dev.off()


MEs = net11$MEs
geneTree = net11$dendrograms[[1]]


# Set up trait file (needs fixing!!!!):
annotation <- metadata
rownames(annotation) <- annotation$File_name
dim(annotation)
names(annotation)
allTraits = annotation[, -c(3, 4)];
dim(allTraits)
names(allTraits)
traitRows = match(rownames(datExpr), allTraits$File_name);
datTraits = allTraits[traitRows, -1, drop = FALSE];
rownames(datTraits) = allTraits[traitRows, 1];
datTraits[,2] <- datTraits[,1]
datTraits[,3] <- datTraits[,1]
datTraits[,4] <- datTraits[,1]
datTraits[,5] <- datTraits[,1]
datTraits[,6] <- datTraits[,1]
names(datTraits) <- c("Day0", "Day1", "Day2", "Day3", "Day4", "Day10")
datTraits[datTraits$Day0==0,1] <- "Day0"
datTraits[datTraits$Day0!="Day0",1] <- 0
datTraits[datTraits$Day0=="Day0",1] <- 1
datTraits[datTraits$Day1==1,2] <- "Day1"
datTraits[datTraits$Day1!="Day1",2] <- 0
datTraits[datTraits$Day1=="Day1",2] <- 1
datTraits[datTraits$Day2==2,3] <- "Day2"
datTraits[datTraits$Day2!="Day2",3] <- 0
datTraits[datTraits$Day2=="Day2",3] <- 1
datTraits[datTraits$Day3==3,4] <- "Day3"
datTraits[datTraits$Day3!="Day3",4] <- 0
datTraits[datTraits$Day3=="Day3",4] <- 1
datTraits[datTraits$Day4==4,5] <- "Day4"
datTraits[datTraits$Day4!="Day4",5] <- 0
datTraits[datTraits$Day4=="Day4",5] <- 1
datTraits[datTraits$Day10==10,6] <- "Day10"
datTraits[datTraits$Day10!="Day10",6] <- 0
datTraits[datTraits$Day10=="Day10",6] <- 1
names(datTraits)



# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p", drop=FALSE)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)


# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
#dim(textMatrix) = dim(moduleTraitCor)



# Display the correlation values within a heatmap plot
plot_path <- "plots/WGCNA_modules_time_p.png"

png(plot_path, height = 600, width = 1000, pointsize = 18)
par(mar = c(6, 10, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


annotation <- metadata
rownames(annotation) <- annotation$File_name
dim(annotation)
names(annotation)
allTraits = annotation[, -c(3, 4)];
dim(allTraits)
names(allTraits)
traitRows = match(rownames(datExpr), allTraits$File_name);
datTraits_time = allTraits[traitRows, -1, drop = FALSE];
rownames(datTraits_time) = allTraits[traitRows, 1];


modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, datTraits_time, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(datTraits_time), sep="");
names(GSPvalue) = paste("p.GS.", names(datTraits_time), sep="");

module = "tan"
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes,]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for timepoint",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
  
colnames(datExpr)[moduleColors=="purple"]
symbol2entrez = match(chicken_dgelist_filt_norm$genes$gene_name, colnames(datExpr))
  
# Create the starting data frame
geneInfo0 <- data.frame(geneSymbol = colnames(datExpr),
                          EntrezID = chicken_dgelist_filt_norm$genes$entrez_gene_id[],
                          moduleColor = moduleColors,
                          geneTraitSignificance,
                          GSPvalue)
# Order modules by their significance for time_point
modOrder <- order(-abs(cor(MEs, datTraits_time, use = "p")));
  
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Timepoint));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "tables/geneInfo_fusion.csv", row.names = FALSE)

geneInfo <- read.delim("tables/geneInfo_fusion.csv", sep = ",", stringsAsFactors = FALSE)
gene_info_mod <- geneInfo
chicken_product_df <- read.delim("chicken_gene_products.tsv", 
                                 header = FALSE, col.names = c("gene_name", "product"))
chicken_product_df <- chicken_product_df[!duplicated(chicken_product_df[,1]),]

egGENENAME <- toTable(org.Gg.egGENENAME)
gene_info_mod$Gene_name <- egGENENAME[match(gene_info_mod$EntrezID, egGENENAME$gene_id),]$gene_name
gene_info_mod$Product <- chicken_product_df[match(gene_info_mod$geneSymbol, chicken_product_df$gene_name),]$product
gene_info_mod$Differentially_expressed <- gene_info_mod$geneSymbol %in% logfc_filt_de_genes_chicken$gene_name
m <- match(gene_info_mod$geneSymbol, rownames(logfc_filt_de_genes_chicken))
gene_info_mod[,43:48] <- logfc_filt_de_genes_chicken[m,]

gene_info_mod <- gene_info_mod[,c(1,2,3,40:48,4:39)]

write.csv(gene_info_mod, file = "tables/geneInfo_fusion_annotated.csv", row.names = FALSE)


# GO and KEGG analyses of the modules

GOenr = GOenrichmentAnalysis(moduleColors, gene_info_mod$entrez_gene_id, organism = "chicken", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
write.table(tab, file = "tables/GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)


# Chicken gene to GO
# GO annotation file from the Gene Ontology Consortium
path_to_chicken_GO_annotation_file <- "goa_chicken.gaf.gz"
chicken_go_annotation <- read.delim(path_to_chicken_GO_annotation_file, stringsAsFactors=FALSE, skip = 31, header = FALSE,
                                    col.names = c("DB", "DB_ID", "Symbol", "Qualifier", "GO_ID", "DB:Reference",
                                                  "Evidence_code", "With_or_from", "Aspect", "DB_Object_Name",
                                                  "DB_Object_Synonyms", "DB_Object_type", "Taxon", "Date", "Assigned_by",
                                                  "Annotation_extension", "Gene_product_form_ID"))
chicken_go <- chicken_go_annotation[,c(3,5)]
go_terms <- toTable(GOTERM)
m <- match(chicken_go$GO_ID, go_terms$go_id)
chicken_go$Term <- go_terms[m,"Term"]
module_go_cats <- chicken_go

module_cat_enrichment <- function(geneInfo, annotation, cat_id, num_genes) {
  # Function for finding the GO enrichment of a module from a WGCNA analysis of chicken and E. tenella data
  # geneInfo should contain gene symbols and module membership
  # annotation should contain gene symbols and  either GO category or KEGG pathway membershp for both chicken and E. tenella genes
  # cat_id is a string that tells what type of categories are being analysed
  # num_genes is the total number of genes in the analysis
  categories <- unique(annotation[,2])
  cat_analysis_results <- as.data.frame(matrix(0, ncol = 5, nrow = length(categories)))
  colnames(cat_analysis_results) <- c(paste(cat_id, "_", geneInfo[1,2], sep = ""), "Term", "N", "num_in_cat", "P_val")
  i <- 0
  while (i < length(categories)) {
    i <- i + 1
    genes_in_cat <- annotation[annotation[,2] == categories[i],]
    N <- dim(genes_in_cat)[1]
    m <- match(genes_in_cat[,1], geneInfo[,1])
    x <- length(m[complete.cases(m)])
    n <- num_genes - N
    k <- dim(geneInfo)[1]
    fisher_result <- fisher.test(matrix(c(x, k-x, N-x, n-(k-x)),nrow=2,ncol=2),alternative="greater")
    cat_results <- data.frame(categories[i], genes_in_cat[1,3], N, x, fisher_result$p.value, stringsAsFactors = FALSE)
    cat_analysis_results[i,] <- cat_results[1,]
  }
  cat_analysis_results <- cat_analysis_results[order(cat_analysis_results$P_val),]
  return(cat_analysis_results)
}

geneInfo_mod <- list()
modNames <- unique(geneInfo$moduleColor)
i <- 0
while (i < length(modNames)) {
  i <- i + 1
  geneInfo_mod[[i]] <- geneInfo[geneInfo$moduleColor == modNames[i],c(1,3)]
}

module_go_cats <- lapply(geneInfo_mod, function(x) module_cat_enrichment(x, chicken_go, "GO_ID", dim(geneInfo)[1]))

i <- 0
while (i < length(module_go_cats)) {
  i <- i + 1
  go_file_name <- paste("tables/wgcna_module_go_kegg/", modNames[i], "_module_go_cats.csv", sep = "")
  write.csv(module_go_cats[[i]], file = go_file_name, row.names = FALSE)
}


# KEGG:
# Chicken gene to KEGG from KEGGrest
chicken_pathway_genes <- keggLink("pathway", "gga")
pathways <- keggList("pathway", "gga")

m <- match(substr(names(chicken_pathway_genes), 5, nchar(names(chicken_pathway_genes))), chicken_dgelist_filt_norm$genes$entrez_gene_id)


chicken_kegg_annotation <- data.frame(gene_symbol = chicken_dgelist$genes$gene_name[m],
                                      kegg_pathway = names(chicken_pathway_genes),
                                      kegg_term = pathways[match(substr(chicken_pathway_genes, 9, nchar(chicken_pathway_genes)),
                                                                 substr(names(pathways), 4, nchar(names(pathways))))],
                                      stringsAsFactors=FALSE)

chicken_kegg_annotation <- chicken_kegg_annotation[complete.cases(chicken_kegg_annotation),]
kegg_annotation <- chicken_kegg_annotation

module_kegg_cats <- lapply(geneInfo_mod, function(x) module_cat_enrichment(x, kegg_annotation, "KEGG_ID", dim(geneInfo)[1]))

i <- 0
while (i < length(module_kegg_cats)) {
  i <- i + 1
  kegg_file_name <- paste("tables/wgcna_module_go_kegg/", modNames[i], "_module_kegg_cats.csv", sep = "")
  write.csv(module_kegg_cats[[i]], file = kegg_file_name, row.names = FALSE)
}

