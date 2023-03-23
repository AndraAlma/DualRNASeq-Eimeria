#To load in all the correct stuff:

install.packages("BiocManager")
BiocManager::install("ggbreak")
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
library("ggbreak")
install_github("vqv/ggbiplot")


path_to_metadata_file <- "~/exjobb/full_metadata.csv"
path_to_chicken_reads <- "~/exjobb/chicken_counts.csv" 

#read in metadata file
metadata <- read.delim(path_to_metadata_file, sep = ";", check.names=FALSE, stringsAsFactors=FALSE)
chicken_metadata <- metadata 
rownames(chicken_metadata) <- chicken_metadata[,1]

# Read in count file 
chicken_reads <- read.delim(path_to_chicken_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)

# create groupings for timepoints
group_chicken <- factor(chicken_metadata[colnames(chicken_reads)[-c(1,2)],"Timepoint"])

# normalize counts to counts per million
chicken_cpm_raw <- cpm(chicken_reads[,-c(1,2)])
rownames(chicken_cpm_raw) <- chicken_reads[,1]

# Check PCA of normalized counts
pca_chicken <- prcomp(chicken_cpm_raw)
ggbiplot(pca_chicken, var.axes = FALSE) + coord_cartesian(xlim = c(0, 160)) +
  geom_text(aes(label = chicken_reads[,1]), hjust = 0.5, vjust = -0.5, size = 3) +
  labs(title = "PCA of normalized chicken read counts per gene") +
  theme(plot.title = element_text(hjust = 0.5))

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

chicken_cpm_raw <- cpm(chicken_reads[,-c(1,2)])
rownames(chicken_cpm_raw) <- chicken_reads[,1]

# Check PCA of normalized counts AGAIN
pca_chicken <- prcomp(chicken_cpm_raw)
ggbiplot(pca_chicken, var.axes = FALSE) + coord_cartesian(xlim = c(0, 80)) +
  geom_text(aes(label = chicken_reads[,1]), hjust = 0.5, vjust = -0.5, size = 3) +
  labs(title = "PCA of normalized chicken read counts per gene") +
  theme(plot.title = element_text(hjust = 0.5))



#Filter out low-expressed genes
chicken_dgelist <- DGEList(counts = chicken_reads[,3:length(chicken_reads)], genes = chicken_reads[,c(1,2)], group = group_chicken)
rownames(chicken_dgelist$counts) <- rownames(chicken_dgelist$genes) <- chicken_dgelist$genes[,1]
keep_chicken <- filterByExpr(chicken_dgelist)
chicken_dgelist_filt <- chicken_dgelist[keep_chicken, ,] 
chicken_dgelist_filt$samples$lib.size <- colSums(chicken_dgelist_filt$counts)


# Normalize gene counts to account for a few genes dominating the counts of the samples
chicken_dgelist_filt_norm <- calcNormFactors(chicken_dgelist_filt)
chicken_labels <- paste(rownames(chicken_dgelist_filt_norm$samples), chicken_dgelist_filt_norm$samples$group, sep = "_")

mds_path <- "~/Exjobb/plots/mds_plots.png"
png(mds_path, height = 1200, width = 800, pointsize = 16)
par(mfrow=c(2,1))
p_chicken <- plotMDS(chicken_dgelist_filt_norm, 
                     labels = chicken_dgelist_filt_norm$samples$group)
dev.off()

chicken_mds_df <- data.frame(x = p_chicken$x, 
                             y = p_chicken$y, 
                             labels = as.character(chicken_dgelist_filt_norm$samples$group), 
                             stringsAsFactors = FALSE)
#i <- 0
#while (i < dim(chicken_mds_df)[1]) {
 # i <- i + 1
#  chicken_mds_df$labels[i] <- substr(chicken_mds_df$labels[i], 3, nchar(chicken_mds_df$labels[i]))
#}

p_chicken_gg <- ggplot(chicken_mds_df, aes(x = x, y = y, color=labels)) +
  geom_point(size = 5) +
  scale_fill_brewer(palette="Dark2")+
  coord_cartesian(ylim = c(-1.8, 1.3)) +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.position = "right")


mds_path <- "~/Exjobb/plots/mds_plots.pdf"
pdf(mds_path, height = 12, width = 12)
p_chicken_gg
dev.off()

chicken_cpm <- cpm.DGEList(chicken_dgelist_filt_norm)
chicken_cpm_df <- data.frame(chicken_cpm, check.names = FALSE)
m <- match(rownames(chicken_cpm_df), chicken_reads$gene_name)
chicken_cpm_df$entrez_gene_id <- chicken_reads$entrez_gene_id[m]
chicken_cpm_df <- chicken_cpm_df[,c(67,1:66)]

write.csv(chicken_cpm_df, file = "~/Exjobb/tables/chicken_counts_norm.csv")

pca_cpm_chicken <- prcomp(t(chicken_cpm))
pca_path <- "~/Exjobb/plots/chicken_pca.png"
png(pca_path, height = 500, width = 520)
ggbiplot(pca_cpm_chicken, var.axes = FALSE, choices = 1:2, alpha = 1) +
  geom_text(aes(label = chicken_dgelist_filt_norm$samples$group), hjust = 0.5, vjust = -0.5, size = 4.5) +
  ggtitle("PCA of normalized chicken read counts per sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()

chicken_dgelist_filt_norm$samples$group <- relevel(chicken_dgelist_filt_norm$samples$group, ref="0")
design_chicken <- model.matrix(~0+group, data = chicken_dgelist_filt_norm$samples)
rownames(design_chicken) <- colnames(chicken_dgelist_filt_norm)


chicken_disp <- estimateDisp(chicken_dgelist_filt_norm, design_chicken, robust = TRUE)
plotBCV(chicken_disp)

fit_chicken <- glmQLFit(chicken_disp, design_chicken, prior.count = 0.125)
#plotQLDisp(fit_chicken)

contrasts_chicken <- makeContrasts(Uvs1=group0-group1, Uvs2=group0-group2, Uvs3=group0-group3, Uvs4=group0-group4, Uvs10=group0-group10,
                                   levels = design_chicken)

#qlf_chicken_B_vs_1 <- glmQLFTest(fit_chicken, coef = 2:7)
#qlf_chicken_U_vs_U.0 <- glmQLFTest(fit_chicken, coef = 8:13)
qlf_chicken_Uvs1 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs1"])
qlf_chicken_Uvs2 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs2"])
qlf_chicken_Uvs3 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs3"])
qlf_chicken_Uvs4 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs4"])
qlf_chicken_Uvs10 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"Uvs10"])

chicken_qlf_list <- list(qlf_chicken_Uvs1,qlf_chicken_Uvs2,qlf_chicken_Uvs3,qlf_chicken_Uvs4,
                         qlf_chicken_Uvs10)
topgenes_chicken_list <- lapply(chicken_qlf_list,
                                function(x) topTags(x, n = dim(chicken_qlf_list[[1]]$table)[1]))

# Define FDR threshold for significance
fdr_threshold <- 0.05
logfc_threshold <- 1

# Volcano plots for all comparisons
timepoints <- c("1_day","2_days","3_days","4_days","10_days")
plot_list <- list()
i = 1
while (i <= length(topgenes_chicken_list)) {
  volcano_path <- paste("~/exjobb/plots/chicken_volcano_lab_", timepoints[i], ".pdf", sep = "")
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
volcano_path <- "~/Exjobb/plots/all_chicken_plots_lab.png"
png(volcano_path, height = 1200, width = 1200)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

volcano_path <- "~/Exjobb/plots/all_chicken_plots_lab.pdf"
pdf(volcano_path, height = 20, width = 20)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()


#Volcano without labels
i = 1
while (i <= length(topgenes_chicken_list)) {
  volcano_path <- paste("~/exjobb/plots/chicken_volcano_", timepoints[i], ".pdf", sep = "")
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
volcano_path <- "~/Exjobb/plots/all_chicken_plots.png"
png(volcano_path, height = 1200, width = 1200)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

volcano_path <- "~/Exjobb/plots/all_chicken_plots.pdf"
pdf(volcano_path, height = 20, width = 20)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

# Compute the number of unique genes differentially expressed at each time point for both datasets, which genes they are and
# which genes are differentially expressed at all time points
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

chicken_disp <- estimateDisp(chicken_dgelist_filt_norm, design_chicken, robust = TRUE)
#eimeria_disp <- estimateDisp(eimeria_dgelist_filt_norm, design_eimeria, robust = TRUE)

plotBCV(chicken_disp)
#plotBCV(eimeria_disp)

# Get a list of all genes that are significantly DE at any timepoint and a list of those that are DE and have a logFC > 1
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

write.table(all_de_genes_chicken, "~/Exjobb/tables/all_de_genes_chicken.csv", sep = ",", row.names = FALSE)

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
logfc_filt_de_genes_chicken$description <- egGENENAME[match(logfc_filt_de_genes_chicken$entrez_gene_id, egGENENAME$gene_id),]$gene_name

write.table(logfc_filt_de_genes_chicken, "~/Exjobb/tables/logfc_filt_de_genes_chicken.csv", sep = ",", row.names = FALSE)

# Get the top 30 most significantly DE genes from each timepoint and export to a table
top_de_genes_path <- "~/Exjobb/tables"
top_de_gene_file_name <- "top_de_genes"
logfc_filtered_topgenes_chicken_list <- lapply(topgenes_chicken_list, function(x) x[abs(x$table$logFC) > logfc_threshold, ])

i <- 0
while (i < length(timepoints)) {
  i <- i + 1
  top_de_genes_file_path_chicken <- paste(top_de_genes_path, "/", top_de_gene_file_name, "_chicken_", timepoints[i], ".csv", sep = "")
  write.table(logfc_filtered_topgenes_chicken_list[[i]][1:30,], top_de_genes_file_path_chicken, sep = ",", row.names = FALSE)
}

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

chicken_fdr_df <- create_gene_heatmap_df(topgenes_chicken_list, "FDR", rownames(topgenes_chicken_list[[1]]$table), timepoints)
chicken_logfc_df <- create_gene_heatmap_df(topgenes_chicken_list, "logFC", rownames(topgenes_chicken_list[[1]]$table), timepoints)

chicken_fdr_df_filt <- filter_gene_heatmap_df(chicken_fdr_df, rownames(logfc_filt_de_genes_chicken))
chicken_logfc_df_filt <- filter_gene_heatmap_df(chicken_logfc_df, rownames(logfc_filt_de_genes_chicken))

heatmap_path <- "~/Exjobb/plots/chicken_genes_pval.pdf"
pdf(heatmap_path, width = 30, height = 240, pointsize = 60)
aspectHeatmap(as.matrix(chicken_fdr_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene p-values", 
              col = colorRampPalette(brewer.pal(9, "YlOrRd"))(500),
              hExp = 1, wExp = 1)
dev.off()

heatmap_path <- "~/Exjobb/plots/chicken_genes_pval.png"
png(heatmap_path, height = 1200, width = 1200)
aspectHeatmap(as.matrix(chicken_fdr_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene p-values", 
              col = colorRampPalette(brewer.pal(9, "YlOrRd"))(500),
              hExp = 1, wExp = 1)
dev.off()

heatmap_path <- "~/Exjobb/plots/chicken_genes_logfc.pdf"
pdf(heatmap_path, width = 240, height = 180, pointsize = 100)
aspectHeatmap(as.matrix(chicken_logfc_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene log fold change", 
              col = colorRampPalette(brewer.pal(9, "RdYlGn"))(500),
              hExp = 1, wExp = 1)
dev.off()

heatmap_path <- "~/Exjobb/plots/chicken_genes_logfc.png"
png(heatmap_path, height = 1200, width = 1200)
aspectHeatmap(as.matrix(chicken_logfc_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene log fold change", 
              col = colorRampPalette(brewer.pal(9, "RdYlGn"))(500),
              hExp = 1, wExp = 1)
dev.off()

# GO category and KEGG pathway analysis of the data.  For each list, the data is in order of increasing time
go_chicken_list <- lapply(chicken_qlf_list, 
                          function(x) goana(x, geneid = chicken_dgelist_filt_norm$genes$entrez_gene_id, FDR = fdr_threshold, species = "Gg"))
#GK <- getGeneKEGGLinks(species.KEGG = "gga")

#keg <- kegga(chicken_qlf_list, species.KEGG="gga", gene.pathway=GK)
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
  go_path <- paste("~/Exjobb/tables/top_go_", timepoints[i], ".csv", sep = "")
  kegg_path <- paste("~/Exjobb/tables/top_kegg_", timepoints[i], ".csv", sep = "")
  write.csv(topgo_chicken_all_list[[i]][1:50,], go_path)
  write.csv(topkegg_chicken_all_list[[i]][1:50,], kegg_path)
}


# GO and KEGG category heatmaps
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

go_chicken_pval_list <- create_cat_pval_df(go_chicken_list, go_chicken_list[[1]]$Term, timepoints)
kegg_chicken_pval_list <- create_cat_pval_df(kegg_chicken_list, kegg_chicken_list[[1]]$Pathway, timepoints)

go_chicken_pval_list_filtered <- lapply(go_chicken_pval_list, function(x) remove_high_pval_rows(x, fdr_threshold/10))
kegg_chicken_pval_list_filtered <- lapply(kegg_chicken_pval_list, function(x) remove_high_pval_rows(x, fdr_threshold))

heatmap_path <- "~/Exjobb/plots/GO_cat_upregulated.pdf"
pdf(heatmap_path, width = 60, height = 180, pointsize = 60)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated GO categories", 
              hExp = 5, wExp = 0.8, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "~/Exjobb/plots/GO_cat_downregulated.pdf"
pdf(heatmap_path, width = 140, height = 180, pointsize = 100)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[2]]), Colv = NA,
              xlab = "Timepoints", main = "Downregulated GO categories", 
              hExp = 2, wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-3,length = 15), seq(-2.95, 2.95, length = 119), seq(3,10,length = 15))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 148)

heatmap_path <- "~/Exjobb/plots/GO_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 60)
heatmap.2(as.matrix(go_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
          xlab = "Timepoints", main = "GO categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45)
dev.off()

heatmap_path <- "~/Exjobb/plots/KEGG_cat_upregulated.pdf"
pdf(heatmap_path, width = 130, height = 120, pointsize = 120)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated KEGG categories", 
              wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "~/Exjobb/plots/KEGG_cat_downregulated.pdf"
pdf(heatmap_path, width = 100, height = 100, pointsize = 100)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[2]]), Colv = NA,
              xlab = "Timepoints", main = "Downregulated KEGG categories", 
              hExp = 2, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-5,length = 11), seq(-4.95, 4.95, length = 199), seq(5,10,length = 11))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 220)

heatmap_path <- "~/Exjobb/plots/KEGG_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 100)
heatmap.2(as.matrix(kegg_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
          xlab = "Timepoints", main = "KEGG categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45)
dev.off()


