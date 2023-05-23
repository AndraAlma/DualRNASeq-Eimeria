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
library("VennDiagram")
library("ggvenn")
library("ggbreak")
install_github("vqv/ggbiplot")


# File paths to count & metadata:
path_to_metadata_file <- "full_metadata_old_new.csv"
path_to_chicken_reads <- "chicken_counts_old_new.csv" 


# Read in files:
metadata <- read.delim(path_to_metadata_file, sep = ";", check.names=FALSE, stringsAsFactors=FALSE)
chicken_metadata <- metadata 
rownames(chicken_metadata) <- chicken_metadata[,1]
chicken_reads <- read.delim(path_to_chicken_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)


# Create groupings for timepoints:
group_chicken <- factor(paste(substr(chicken_metadata[colnames(chicken_reads)[-c(1,2)],"Infection_status"],1,1),
                              chicken_metadata[colnames(chicken_reads)[-c(1,2)],"Timepoint"],
                              sep="."))


# Normalize counts to counts per million:
chicken_cpm_raw <- cpm(chicken_reads[,-c(1,2)])
rownames(chicken_cpm_raw) <- chicken_reads[,1]


# Check PCA of normalized counts to detect outliers:
pca_chicken <- prcomp(chicken_cpm_raw)
pca_path <- "plots/old_and_new/pca.png"
png(pca_path, height = 1200, width = 1200)
ggbiplot(pca_chicken, var.axes = FALSE) + coord_cartesian(xlim = c(0, 160)) +
  geom_text(aes(label = chicken_reads[,1]), hjust = 0.5, vjust = -0.5, size = 3) +
  labs(title = "PCA of normalized chicken read counts per gene") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


# Remove outliers:
chicken_remove = which(chicken_reads[,1] == "LOC112533599")
chicken_reads = chicken_reads[-chicken_remove,]

chicken_remove = which(chicken_reads[,1] == "LOC112533601") 
chicken_reads = chicken_reads[-chicken_remove,]

chicken_remove = which(chicken_reads[,1] == "EYY68_mgr01")
chicken_reads = chicken_reads[-chicken_remove,]

chicken_remove = which(chicken_reads[,1] == "EYY68_mgr02")
chicken_reads = chicken_reads[-chicken_remove,]


# Check PCA of normalized counts again:
chicken_cpm_raw <- cpm(chicken_reads[,-c(1,2)])
rownames(chicken_cpm_raw) <- chicken_reads[,1]

pca_chicken <- prcomp(chicken_cpm_raw)
pca_path <- "plots/old_and_new/pca_checked.png"
png(pca_path, height = 1200, width = 1200)
ggbiplot(pca_chicken, var.axes = FALSE) + coord_cartesian(xlim = c(0, 80)) +
  geom_text(aes(label = chicken_reads[,1]), hjust = 0.5, vjust = -0.5, size = 3) +
  labs(title = "PCA of normalized chicken read counts per gene") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


# Filter out low-expressed genes:
chicken_dgelist <- DGEList(counts = chicken_reads[,3:length(chicken_reads)], genes = chicken_reads[,c(1,2)], group = group_chicken)
rownames(chicken_dgelist$counts) <- rownames(chicken_dgelist$genes) <- chicken_dgelist$genes[,1]
keep_chicken <- filterByExpr(chicken_dgelist)
chicken_dgelist_filt <- chicken_dgelist[keep_chicken, ,] 
chicken_dgelist_filt$samples$lib.size <- colSums(chicken_dgelist_filt$counts)


# Normalize gene counts to account for a few genes dominating the counts of the samples:
chicken_dgelist_filt_norm <- calcNormFactors(chicken_dgelist_filt)
chicken_labels <- paste(rownames(chicken_dgelist_filt_norm$samples), chicken_dgelist_filt_norm$samples$group, sep = "_")


# Normalize to counts per million, add entrez gene IDs and print normalized counts to table:
chicken_cpm <- cpm.DGEList(chicken_dgelist_filt_norm)
chicken_cpm_df <- data.frame(chicken_cpm, check.names = FALSE)
m <- match(rownames(chicken_cpm_df), chicken_reads$gene_name)
chicken_cpm_df$entrez_gene_id <- chicken_reads$entrez_gene_id[m]
chicken_cpm_df <- chicken_cpm_df[,c(103,1:102)]

write.csv(chicken_cpm_df, file = "tables/old_and_new/chicken_counts_norm.csv")


# Plot PCA:
pca_cpm_chicken <- prcomp(t(chicken_cpm))
pca_p <- pca_cpm_chicken
label_p <- chicken_dgelist_filt_norm$samples$group
pca_path <- "plots/old_and_new/chicken_pca.png"
png(pca_path, height = 1000, width = 1000)
ggbiplot(pca_cpm_chicken, var.axes = FALSE, choices = 1:2, alpha = 1) +
  geom_text(aes(label = chicken_dgelist_filt_norm$samples$group), hjust = 0.5, vjust = -0.5, size = 4.5) +
  ggtitle("PCA of normalized chicken read counts per sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()



# Comparing primary to secondary dataset:


# Create design matrix:
chicken_dgelist_filt_norm$samples$group <- relevel(chicken_dgelist_filt_norm$samples$group, ref="p.0")
design_chicken <- model.matrix(~0+group, data = chicken_dgelist_filt_norm$samples)
rownames(design_chicken) <- colnames(chicken_dgelist_filt_norm)

chicken_disp <- estimateDisp(chicken_dgelist_filt_norm, design_chicken, robust = TRUE)
plotBCV(chicken_disp)

# Create contrasts between primary & secondary infections for each day:
fit_chicken <- glmQLFit(chicken_disp, design_chicken, prior.count = 0.125)
contrasts_chicken <- makeContrasts(PS.1=groups.1-groupp.1, PS.2=groups.2-groupp.2,PS.3=groups.3-groupp.3,PS.4=groups.4-groupp.4,PS.10=groups.10-groupp.10,
                                   levels = design_chicken)

qlf_chicken_PS1 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"PS.1"])
qlf_chicken_PS2 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"PS.2"])
qlf_chicken_PS3 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"PS.3"])
qlf_chicken_PS4 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"PS.4"])
qlf_chicken_PS10 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"PS.10"])

chicken_qlf_list <- list(qlf_chicken_PS1,qlf_chicken_PS2,qlf_chicken_PS3,qlf_chicken_PS4,
                         qlf_chicken_PS10)
topgenes_chicken_list <- lapply(chicken_qlf_list,
                                function(x) topTags(x, n = dim(chicken_qlf_list[[1]]$table)[1]))


# Define FDR & logFC threshold for significance:
fdr_threshold <- 0.05
logfc_threshold <- 1


# Create volcano plots for all comparisons:
timepoints <- c("1_dpi","2_dpi","3_dpi","4_dpi","10_dpi")
plot_list <- list()
i = 1
while (i <= length(topgenes_chicken_list)) {
  volcano_path <- paste("plots/old_and_new/chicken_volcano_lab_", timepoints[i], ".pdf", sep = "")
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
volcano_path <- "plots/old_and_new/all_chicken_plots_lab.png"
png(volcano_path, height = 1200, width = 1200)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

volcano_path <- "plots/old_and_new/all_chicken_plots_lab.pdf"
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

write.table(all_de_genes_chicken, "tables/old_and_new/all_de_genes_chicken.csv", sep = ",", row.names = FALSE)

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
write.table(logfc_filt_de_genes_chicken, "tables/old_and_new/logfc_filt_de_genes_chicken.csv", sep = ",", row.names = FALSE)


# Get the top 30 most significantly DE genes from each timepoint and export to a table:
top_de_genes_path <- "tables/old_and_new"
top_de_gene_file_name <- "top_de_genes"
logfc_filtered_topgenes_chicken_list <- lapply(topgenes_chicken_list, function(x) x[abs(x$table$logFC) > logfc_threshold, ])

i <- 0
while (i < length(timepoints)) {
  i <- i + 1
  top_de_genes_file_path_chicken <- paste(top_de_genes_path, "/", top_de_gene_file_name, "_chicken_", timepoints[i], ".csv", sep = "")
  write.table(logfc_filtered_topgenes_chicken_list[[i]][1:30,], top_de_genes_file_path_chicken, sep = ",", row.names = FALSE)
}


# Filter for FDR as well, and print top 40 to a list:
top_de_genes_path <- "tables/old_and_new"
top_de_gene_file_name <- "top_de_genes_FDR"
i <- 0
while (i < length(timepoints)) {
  i <- i + 1
  top_de_genes_file_path_chicken <- paste(top_de_genes_path, "/", top_de_gene_file_name, "_chicken_", timepoints[i], ".csv", sep = "")
  write.table(de_genes_chicken_all[[i]][1:40,], top_de_genes_file_path_chicken, sep = ",", row.names = FALSE)
}


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
  go_path <- paste("tables/old_and_new/top_go_", timepoints[i], ".csv", sep = "")
  kegg_path <- paste("tables/old_and_new/top_kegg_", timepoints[i], ".csv", sep = "")
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

heatmap_path <- "plots/old_and_new/GO_cat_upregulated.pdf"
pdf(heatmap_path, width = 60, height = 180, pointsize = 60)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated GO categories", 
              hExp = 5, wExp = 0.8, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "plots/old_and_new/GO_cat_downregulated.pdf"
pdf(heatmap_path, width = 140, height = 180, pointsize = 100)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[2]]), Colv = NA,
              xlab = "Timepoints", main = "Downregulated GO categories", 
              hExp = 2, wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-3,length = 15), seq(-2.95, 2.95, length = 119), seq(3,10,length = 15))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 148)

heatmap_path <- "plots/old_and_new/GO_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 60)
heatmap.2(as.matrix(go_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
          xlab = "Timepoints", main = "GO categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45)
dev.off()

heatmap_path <- "plots/old_and_new/KEGG_cat_upregulated.pdf"
pdf(heatmap_path, width = 130, height = 120, pointsize = 120)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated KEGG categories", 
              wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "plots/old_and_new/KEGG_cat_downregulated.pdf"
pdf(heatmap_path, width = 100, height = 100, pointsize = 100)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[2]]), Colv = NA,
              xlab = "Timepoints", main = "Downregulated KEGG categories", 
              hExp = 2, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-5,length = 11), seq(-4.95, 4.95, length = 199), seq(5,10,length = 11))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 220)

heatmap_path <- "plots/old_and_new/KEGG_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 100)
heatmap.2(as.matrix(kegg_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
          xlab = "Timepoints", main = "KEGG categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45)
dev.off()


# Investigate immune gene groups:
cyto_chemo <- c("CCL19", "CCL4-2", "CCL26", "IL8L2", "IL8L1", "CNMD")
IFN_down <- c("CDH23", "LOC771422", "CLCF1", "CYP1A1", "LOC112532392", "DNAH1",
              "LOC112532392", "GPNMB")
IFN_up <- c("ABCB1LA", "LYGL", "LOC107054696", "GVINP1", "GBP", "VSIG1", 
            "LOC100858381", "BATF3", "IFI6", "IFIT5", "MX1", "RSAD2", "DDX4", 
            "AMY2A", "USP18", "LOC100858381", "CMPK2", "LY6E", "HBA1", "KCNG3", 
            "LOC112532459", "VSIG1", "SGCZ", "SLC27A6")
# Split up in 2 smaller subsets:
IFN_up1 <- c("LYGL", "LY6E", "IFI6", "IFIT5", "RSAD2", "USP18", "MX1", "CMPK2",
             "AMY2A", "DDX4",  "SLC27A6", "VSIG1")
IFN_up2 <- c("GBP", "LOC107054696", "GVINP1", "ABCB1LA","SGCZ",
             "LOC112532459", "LOC100858381", "BATF3", "HBA1")

# Function to create line plots:
plot_de_cat_genes <- function(de_cat_genes, experimental_conditions, plot_title, fdr_thresh = 0.05, 
                              logfc_thresh = 1, num_samples_fdr_thresh = 1, num_samples_logfc_thresh = 1,
                              cat_plot_path, plot_scale = c(-5,4)) {
  # Produces a line plot of log2 fold change across a set of samples for genes in a specific
  # GO or KEGG category
  i <- 1
  j <- 1
  de_cat_gene_df <- as.data.frame(matrix(0, 
                                         ncol = dim(de_cat_genes[[1]])[2], 
                                         nrow = dim(de_cat_genes[[1]])[1]*length(de_cat_genes)))
  colnames(de_cat_gene_df) <- colnames(de_cat_genes[[1]])
  num_conditions <- length(de_cat_genes)
  num_genes <- dim(de_cat_genes[[1]])[1]
  while (i <= num_genes) {
    k <- 1
    while (k <= num_conditions) {
      de_cat_gene_df[j,] <- de_cat_genes[[k]][i,]
      j <- j + 1
      k <- k + 1
    }
    i <- i + 1
  }
  
  below_fdr_threshold <- de_cat_gene_df$FDR < fdr_thresh
  de_cat_gene_df$below_fdr_threshold <- below_fdr_threshold
  curr_num_genes <- dim(de_cat_genes[[1]])[1]
  neg_logfc_thresh <- -logfc_thresh
  
  i <- 0
  while (i < curr_num_genes) {
    gene_start <- num_conditions*i + 1
    gene_end <- num_conditions*(i + 1)
    gene_below_thresh <- FALSE
    j <- num_conditions*i + 1
    num_samples_below_fdr_thresh = 0
    num_samples_above_logfc_thresh = 0
    while (j <= num_conditions*(i+1)) {
      if (de_cat_gene_df$below_fdr_threshold[j]) {
        num_samples_below_fdr_thresh <- num_samples_below_fdr_thresh + 1
      }
      if (de_cat_gene_df$logFC[j] >= logfc_thresh || de_cat_gene_df$logFC[j] <= neg_logfc_thresh) {
        num_samples_above_logfc_thresh <- num_samples_above_logfc_thresh + 1
      }
      j <- j + 1
    }
    if (num_samples_below_fdr_thresh >= num_samples_fdr_thresh && num_samples_above_logfc_thresh >= num_samples_logfc_thresh) {
      i <- i + 1
    } 
    else {
      gene_start <- num_conditions*i + 1
      gene_end <- num_conditions*(i + 1)
      de_cat_gene_df <- de_cat_gene_df[-c(gene_start:gene_end),]
    }
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
  }
  
  if (curr_num_genes > 20) {
    # To increase readability, the plot is split into two if there are too many genes being plotted
    split_vec_1 <- rep(1, curr_num_genes/2)
    split_vec_2 <- rep(2, curr_num_genes/2)
    split_vec <- c(split_vec_1, split_vec_1, split_vec_1, split_vec_1, split_vec_1, 
                   split_vec_2, split_vec_2, split_vec_2, split_vec_2, split_vec_2)
    if (curr_num_genes*5 > length(split_vec)) {
      split_vec <- c(split_vec, rep(2, num_conditions))
    }
    de_cat_gene_df_list <- split(de_cat_gene_df, split_vec)
    de_cat_gene_df <- de_cat_gene_df_list[[1]]
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
    
    de_cat_gene_df$gene_timepoints <- rep(experimental_conditions, curr_num_genes)
    cat_plot_path <- gsub(".png", "_1.png", cat_plot_path)
    png(cat_plot_path, width = 1200, height = 800)
    p <- ggplot(de_cat_gene_df, aes(x=gene_timepoints, y=logFC, group=gene_name)) +
      geom_line(aes(color = gene_name), size = 1.25) +
      geom_point(aes(shape=below_fdr_threshold, color = gene_name), size = 5) +
      scale_shape_manual(values = c(1,17)) +
      scale_x_discrete(limits=experimental_conditions) +
      theme_bw() +
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5, size = 36),
            legend.title = element_text(size = 28),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 24)) +
      guides(color = guide_legend(order = 1),
             shape = guide_legend(order = 2)) +
      coord_cartesian(ylim = plot_scale, xlim = c(0.25, 10.5)) +
      labs(shape = "Below FDR threshold", color = "Gene symbol") +
      xlab("Time post-infection [days]") + ylab("log2(Fold change)") +
      geom_hline(yintercept = 0)
    print(p)
    dev.off()
    cat_plot_path <- gsub("_1.png", "_2.png", cat_plot_path)
    
    de_cat_gene_df <- de_cat_gene_df_list[[2]]
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
  }
  
  de_cat_gene_df$gene_timepoints <- rep(experimental_conditions, curr_num_genes)
  png(cat_plot_path, width = 1200, height = 800)
  par(mar = c(5, 5, 5, 10))
  p <- ggplot(de_cat_gene_df, aes(x=gene_timepoints, y=logFC, color=gene_name)) +
    scale_colour_brewer(palette="Paired")+
    geom_line(aes(color = gene_name), size = 1.25) +
    geom_point(aes(shape=below_fdr_threshold, color = gene_name), size = 5) +
    scale_shape_manual(values = c(1,17)) +
    scale_x_discrete(limits=experimental_conditions) +
    theme_bw() +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 36),
          legend.title = element_text(size = 23),
          legend.text = element_text(size = 18),
          legend.position = "bottom",
          legend.box = "vertical",
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    coord_cartesian(ylim = plot_scale, xlim = c(0.25, 10.5)) +
    labs(shape = "Below FDR threshold", color = "Gene symbol") +
    xlab("Time post-infection [days]") + ylab("log2(Fold change)") +
    geom_hline(yintercept = 0)
  print(p)
  dev.off()
}


# Plot each group:
cat_gene_list <- lapply(topgenes_chicken_list, function(x) x$table[match(cyto_chemo, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list, c(1,2,3,4,10), "Cytokines & Chemokines", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-5,2.5), "plots/old_and_new/cyto_chemo_expr.png")
cat_gene_list <- lapply(topgenes_chicken_list, function(x) x$table[match(IFN_down, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list, c(1,2,3,4,10), "IFN genes, downregulated", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-6,7.5), "plots/old_and_new/IFN_down_expr.png")
cat_gene_list <- lapply(topgenes_chicken_list, function(x) x$table[match(IFN_up, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list, c(1,2,3,4,10), "IFN genes, upregulated", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-4,5.25), "plots/old_and_new/IFN_up_expr.png")
cat_gene_list <- lapply(topgenes_chicken_list, function(x) x$table[match(IFN_up1, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list, c(1,2,3,4,10), "IFN genes, upregulated", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-5,3), "plots/old_and_new/IFN_up_expr_col1.png")
cat_gene_list <- lapply(topgenes_chicken_list, function(x) x$table[match(IFN_up2, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list, c(1,2,3,4,10), "IFN genes, upregulated", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-4,6), "plots/old_and_new/IFN_up_expr_col2.png")




# Only primary infection data:
# Create design matrix:
chicken_dgelist_filt_norm$samples$group <- relevel(chicken_dgelist_filt_norm$samples$group, ref="p.0")
design_chicken <- model.matrix(~0+group, data = chicken_dgelist_filt_norm$samples)
rownames(design_chicken) <- colnames(chicken_dgelist_filt_norm)

chicken_disp <- estimateDisp(chicken_dgelist_filt_norm, design_chicken, robust = TRUE)
plotBCV(chicken_disp)

# Create contrasts between reference timepoint and all other timepoints for DE:
fit_chicken <- glmQLFit(chicken_disp, design_chicken, prior.count = 0.125)
contrasts_chicken <- makeContrasts(D.1=groupp.1-groupp.0, D.2=groupp.2-groupp.0,D.3=groupp.3-groupp.0,D.4=groupp.4-groupp.0,D.10=groupp.10-groupp.0,
                                   levels = design_chicken)
qlf_chicken_D1 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"D.1"])
qlf_chicken_D2 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"D.2"])
qlf_chicken_D3 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"D.3"])
qlf_chicken_D4 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"D.4"])
qlf_chicken_D10 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"D.10"])

chicken_qlf_list <- list(qlf_chicken_D1,qlf_chicken_D2,qlf_chicken_D3,qlf_chicken_D4,
                         qlf_chicken_D10)
topgenes_chicken_list <- lapply(chicken_qlf_list,
                                function(x) topTags(x, n = dim(chicken_qlf_list[[1]]$table)[1]))


# Define FDR & logFC threshold for significance:
fdr_threshold <- 0.05
logfc_threshold <- 1


# Create volcano plots for all comparisons:
timepoints <- c("1_dpi","2_dpi","3_dpi","4_dpi","10_dpi")
plot_list <- list()
i = 1
while (i <= length(topgenes_chicken_list)) {
  volcano_path <- paste("plots/old_and_new/old_chicken_volcano_lab_", timepoints[i], ".pdf", sep = "")
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

volcano_path <- "plots/old_and_new/old_all_chicken_plots_lab.png"
png(volcano_path, height = 1200, width = 1200)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

volcano_path <- "plots/old_and_new/old_all_chicken_plots_lab.pdf"
pdf(volcano_path, height = 20, width = 20)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

# Volcano plots without labels:
i = 1
while (i <= length(topgenes_chicken_list)) {
  volcano_path <- paste("plots/old_and_new/old_chicken_volcano_", timepoints[i], ".pdf", sep = "")
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
    theme(plot.title = element_text(hjust = 0.5, size = 24),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 20))
  plot_list[[i]] <- p
  print(p)
  dev.off()
  i = i + 1
}

# Save plot list for comaprison with secondary infection:
plot_list_primary <- plot_list 

volcano_path <- "plots/old_and_new/old_all_chicken_plots.png"
png(volcano_path, height = 1200, width = 1200)
multiplot(plotlist = plot_list, layout = matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

volcano_path <- "plots/old_and_new/old_all_chicken_plots.pdf"
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

write.table(all_de_genes_chicken, "tables/old_and_new/old_all_de_genes_chicken.csv", sep = ",", row.names = FALSE)

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
write.table(logfc_filt_de_genes_chicken, "tables/old_and_new/old_logfc_filt_de_genes_chicken.csv", sep = ",", row.names = FALSE)


# Get the top 30 most significantly DE genes from each timepoint and export to a table:
top_de_genes_path <- "tables/old_and_new"
top_de_gene_file_name <- "old_top_de_genes"
logfc_filtered_topgenes_chicken_list <- lapply(topgenes_chicken_list, function(x) x[abs(x$table$logFC) > logfc_threshold, ])

i <- 0
while (i < length(timepoints)) {
  i <- i + 1
  top_de_genes_file_path_chicken <- paste(top_de_genes_path, "/", top_de_gene_file_name, "_chicken_", timepoints[i], ".csv", sep = "")
  write.table(logfc_filtered_topgenes_chicken_list[[i]][1:30,], top_de_genes_file_path_chicken, sep = ",", row.names = FALSE)
}


# Filter for FDR as well, and print top 40 to a list:
top_de_genes_path <- "tables/old_and_new"
top_de_gene_file_name <- "old_top_de_genes_FDR"
i <- 0
while (i < length(timepoints)) {
  i <- i + 1
  top_de_genes_file_path_chicken <- paste(top_de_genes_path, "/", top_de_gene_file_name, "_chicken_", timepoints[i], ".csv", sep = "")
  write.table(de_genes_chicken_all[[i]][1:40,], top_de_genes_file_path_chicken, sep = ",", row.names = FALSE)
}


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
  go_path <- paste("tables/old_and_new/old_top_go_", timepoints[i], ".csv", sep = "")
  kegg_path <- paste("tables/old_and_new/old_top_kegg_", timepoints[i], ".csv", sep = "")
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

# Drop days from dataframe without any significant genes:
DE_go_chicken_list <- go_chicken_list[-1]
DE_kegg_chicken_list <- kegg_chicken_list[-1]


# Create heatmaps for GO and KEGG categories:
DE_timepoints <-timepoints[2:length(timepoints)]
go_chicken_pval_list <- create_cat_pval_df(DE_go_chicken_list, DE_go_chicken_list[[1]]$Term, DE_timepoints)
kegg_chicken_pval_list <- create_cat_pval_df(DE_kegg_chicken_list, DE_kegg_chicken_list[[1]]$Pathway, DE_timepoints)

go_chicken_pval_list_filtered <- lapply(go_chicken_pval_list, function(x) remove_high_pval_rows(x, fdr_threshold/10))
kegg_chicken_pval_list_filtered <- lapply(kegg_chicken_pval_list, function(x) remove_high_pval_rows(x, fdr_threshold))

heatmap_path <- "plots/old_and_new/old_GO_cat_upregulated.pdf"
pdf(heatmap_path, width = 60, height = 180, pointsize = 60)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated GO categories", 
              hExp = 5, wExp = 0.8, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "plots/old_and_new/old_GO_cat_downregulated.pdf"
pdf(heatmap_path, width = 140, height = 180, pointsize = 100)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[2]]), Colv = NA,
              xlab = "Timepoints", main = "Downregulated GO categories", 
              hExp = 2, wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-3,length = 15), seq(-2.95, 2.95, length = 119), seq(3,10,length = 15))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 148)

heatmap_path <- "plots/old_and_new/old_GO_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 60)
heatmap.2(as.matrix(go_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
          xlab = "Timepoints", main = "GO categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45)
dev.off()

heatmap_path <- "plots/old_and_new/old_KEGG_cat_upregulated.pdf"
pdf(heatmap_path, width = 130, height = 120, pointsize = 120)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated KEGG categories", 
              wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "plots/old_and_new/old_KEGG_cat_downregulated.pdf"
pdf(heatmap_path, width = 100, height = 100, pointsize = 100)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[2]]), Colv = NA,
              xlab = "Timepoints", main = "Downregulated KEGG categories", 
              hExp = 2, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-5,length = 11), seq(-4.95, 4.95, length = 199), seq(5,10,length = 11))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 220)

heatmap_path <- "plots/old_and_new/old_KEGG_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 100)
heatmap.2(as.matrix(kegg_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
          xlab = "Timepoints", main = "KEGG categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45)
dev.off()


# Investigate immune gene groups:
cyto_chemo <- c("CCL19", "CCL4-2", "CCL26", "IL8L2", "IL8L1", "CNMD")
IFN_down <- c("CDH23", "LOC771422", "CLCF1", "CYP1A1", "LOC112532392", "DNAH1",
              "LOC112532392", "GPNMB")
IFN_up <- c("ABCB1LA", "LYGL", "LOC107054696", "GVINP1", "GBP", "VSIG1", 
            "LOC100858381", "BATF3", "IFI6", "IFIT5", "MX1", "RSAD2", "DDX4", 
            "AMY2A", "USP18", "LOC100858381", "CMPK2", "LY6E", "HBA1", "KCNG3", 
            "LOC112532459", "VSIG1", "SGCZ", "SLC27A6")

# split up into 2 smaller subsets:
IFN_up1 <- c("LYGL", "LY6E", "IFI6", "IFIT5", "RSAD2", "USP18", "MX1", "CMPK2",
             "AMY2A", "DDX4",  "SLC27A6", "VSIG1")
IFN_up2 <- c("GBP", "LOC107054696", "GVINP1", "ABCB1LA","SGCZ",
             "LOC112532459", "LOC100858381", "BATF3", "HBA1")


# Function to create plots:
plot_de_cat_genes <- function(de_cat_genes, experimental_conditions, plot_title, fdr_thresh = 0.05, 
                              logfc_thresh = 1, num_samples_fdr_thresh = 1, num_samples_logfc_thresh = 1,
                              cat_plot_path, plot_scale = c(-5,4)) {
  # Produces a line plot of log2 fold change across a set of samples for genes in a specific
  # GO or KEGG category
  i <- 1
  j <- 1
  de_cat_gene_df <- as.data.frame(matrix(0, 
                                         ncol = dim(de_cat_genes[[1]])[2], 
                                         nrow = dim(de_cat_genes[[1]])[1]*length(de_cat_genes)))
  colnames(de_cat_gene_df) <- colnames(de_cat_genes[[1]])
  num_conditions <- length(de_cat_genes)
  num_genes <- dim(de_cat_genes[[1]])[1]
  while (i <= num_genes) {
    k <- 1
    while (k <= num_conditions) {
      de_cat_gene_df[j,] <- de_cat_genes[[k]][i,]
      j <- j + 1
      k <- k + 1
    }
    i <- i + 1
  }
  
  below_fdr_threshold <- de_cat_gene_df$FDR < fdr_thresh
  de_cat_gene_df$below_fdr_threshold <- below_fdr_threshold
  curr_num_genes <- dim(de_cat_genes[[1]])[1]
  neg_logfc_thresh <- -logfc_thresh
  
  i <- 0
  while (i < curr_num_genes) {
    gene_start <- num_conditions*i + 1
    gene_end <- num_conditions*(i + 1)
    gene_below_thresh <- FALSE
    j <- num_conditions*i + 1
    num_samples_below_fdr_thresh = 0
    num_samples_above_logfc_thresh = 0
    while (j <= num_conditions*(i+1)) {
      if (de_cat_gene_df$below_fdr_threshold[j]) {
        num_samples_below_fdr_thresh <- num_samples_below_fdr_thresh + 1
      }
      if (de_cat_gene_df$logFC[j] >= logfc_thresh || de_cat_gene_df$logFC[j] <= neg_logfc_thresh) {
        num_samples_above_logfc_thresh <- num_samples_above_logfc_thresh + 1
      }
      j <- j + 1
    }
    if (num_samples_below_fdr_thresh >= num_samples_fdr_thresh && num_samples_above_logfc_thresh >= num_samples_logfc_thresh) {
      i <- i + 1
    } 
    else {
      gene_start <- num_conditions*i + 1
      gene_end <- num_conditions*(i + 1)
      de_cat_gene_df <- de_cat_gene_df[-c(gene_start:gene_end),]
    }
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
  }
  
  if (curr_num_genes > 20) {
    # To increase readability, the plot is split into two if there are too many genes being plotted
    split_vec_1 <- rep(1, curr_num_genes/2)
    split_vec_2 <- rep(2, curr_num_genes/2)
    split_vec <- c(split_vec_1, split_vec_1, split_vec_1, split_vec_1, split_vec_1, 
                   split_vec_2, split_vec_2, split_vec_2, split_vec_2, split_vec_2)
    if (curr_num_genes*5 > length(split_vec)) {
      split_vec <- c(split_vec, rep(2, num_conditions))
    }
    de_cat_gene_df_list <- split(de_cat_gene_df, split_vec)
    de_cat_gene_df <- de_cat_gene_df_list[[1]]
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
    
    de_cat_gene_df$gene_timepoints <- rep(experimental_conditions, curr_num_genes)
    cat_plot_path <- gsub(".png", "_1.png", cat_plot_path)
    png(cat_plot_path, width = 1200, height = 800)
    p <- ggplot(de_cat_gene_df, aes(x=gene_timepoints, y=logFC, group=gene_name)) +
      geom_line(aes(color = gene_name), size = 1.25) +
      geom_point(aes(shape=below_fdr_threshold, color = gene_name), size = 5) +
      scale_shape_manual(values = c(1,17)) +
      scale_x_discrete(limits=experimental_conditions) +
      theme_bw() +
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5, size = 36),
            legend.title = element_text(size = 28),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 24)) +
      guides(color = guide_legend(order = 1),
             shape = guide_legend(order = 2)) +
      coord_cartesian(ylim = plot_scale, xlim = c(0.25, 10.5)) +
      labs(shape = "Below FDR threshold", color = "Gene symbol") +
      xlab("Time post-infection [days]") + ylab("log2(Fold change)") +
      geom_hline(yintercept = 0)
    print(p)
    dev.off()
    cat_plot_path <- gsub("_1.png", "_2.png", cat_plot_path)
    
    de_cat_gene_df <- de_cat_gene_df_list[[2]]
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
  }
  
  de_cat_gene_df$gene_timepoints <- rep(experimental_conditions, curr_num_genes)
  png(cat_plot_path, width = 1200, height = 800)
  par(mar = c(5, 5, 5, 10))
  p <- ggplot(de_cat_gene_df, aes(x=gene_timepoints, y=logFC, color=gene_name)) +
    scale_colour_brewer(palette="Paired")+
    geom_line(aes(color = gene_name), size = 1.25) +
    geom_point(aes(shape=below_fdr_threshold, color = gene_name), size = 5) +
    scale_shape_manual(values = c(1,17)) +
    scale_x_discrete(limits=experimental_conditions) +
    theme_bw() +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 36),
          legend.title = element_text(size = 23),
          legend.text = element_text(size = 18),
          legend.position = "bottom",
          legend.box = "vertical",
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    coord_cartesian(ylim = plot_scale, xlim = c(0.25, 10.5)) +
    labs(shape = "Below FDR threshold", color = "Gene symbol") +
    xlab("Time post-infection [days]") + ylab("log2(Fold change)") +
    geom_hline(yintercept = 0)
  print(p)
  dev.off()
}




# Plot each group:
cat_gene_list_p_cc <- lapply(topgenes_chicken_list, function(x) x$table[match(cyto_chemo, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list_p_cc, c(1,2,3,4,10), "Cytokines & Chemokines, primary infection", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-2.5,4), "plots/old_and_new/old_cyto_chemo_expr.png")
cat_gene_list_p_d <- lapply(topgenes_chicken_list, function(x) x$table[match(IFN_down, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list_p_d, c(1,2,3,4,10), "IFN genes downregulated, primary infection ", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-3,7.5), "plots/old_and_new/old_IFN_down_expr.png")
cat_gene_list_p_u <- lapply(topgenes_chicken_list, function(x) x$table[match(IFN_up, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list_p_u, c(1,2,3,4,10), "IFN genes upregulated, primary infection", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-4,5.25), "plots/old_and_new/old_IFN_up_expr.png")
cat_gene_list_p_u1 <- lapply(topgenes_chicken_list, function(x) x$table[match(IFN_up1, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list_p_u1, c(1,2,3,4,10), "IFN genes upregulated, primary infection", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-3,6), "plots/old_and_new/old_IFN_up_expr_col1.png")
cat_gene_list_p_u2 <- lapply(topgenes_chicken_list, function(x) x$table[match(IFN_up2, x$table$gene_name),])
plot_de_cat_genes(cat_gene_list_p_u2, c(1,2,3,4,10), "IFN genes upregulated, primary infection", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-4,6), "plots/old_and_new/old_IFN_up_expr_col2.png")


# Line plots side by side for primary and secondary:
plot_de_cat_genes2 <- function(de_cat_genes, experimental_conditions, plot_title, fdr_thresh = 0.05, 
                              logfc_thresh = 1, num_samples_fdr_thresh = 1, num_samples_logfc_thresh = 1,
                              plot_scale = c(-5,4)) {
  # Produces a line plot of log2 fold change across a set of samples for genes in a specific
  # GO or KEGG category
  i <- 1
  j <- 1
  de_cat_gene_df <- as.data.frame(matrix(0, 
                                         ncol = dim(de_cat_genes[[1]])[2], 
                                         nrow = dim(de_cat_genes[[1]])[1]*length(de_cat_genes)))
  colnames(de_cat_gene_df) <- colnames(de_cat_genes[[1]])
  num_conditions <- length(de_cat_genes)
  num_genes <- dim(de_cat_genes[[1]])[1]
  while (i <= num_genes) {
    k <- 1
    while (k <= num_conditions) {
      de_cat_gene_df[j,] <- de_cat_genes[[k]][i,]
      j <- j + 1
      k <- k + 1
    }
    i <- i + 1
  }
  
  below_fdr_threshold <- de_cat_gene_df$FDR < fdr_thresh
  de_cat_gene_df$below_fdr_threshold <- below_fdr_threshold
  curr_num_genes <- dim(de_cat_genes[[1]])[1]
  neg_logfc_thresh <- -logfc_thresh
  
  i <- 0
  while (i < curr_num_genes) {
    gene_start <- num_conditions*i + 1
    gene_end <- num_conditions*(i + 1)
    gene_below_thresh <- FALSE
    j <- num_conditions*i + 1
    num_samples_below_fdr_thresh = 0
    num_samples_above_logfc_thresh = 0
    while (j <= num_conditions*(i+1)) {
      if (de_cat_gene_df$below_fdr_threshold[j]) {
        num_samples_below_fdr_thresh <- num_samples_below_fdr_thresh + 1
      }
      if (de_cat_gene_df$logFC[j] >= logfc_thresh || de_cat_gene_df$logFC[j] <= neg_logfc_thresh) {
        num_samples_above_logfc_thresh <- num_samples_above_logfc_thresh + 1
      }
      j <- j + 1
    }
    if (num_samples_below_fdr_thresh >= num_samples_fdr_thresh && num_samples_above_logfc_thresh >= num_samples_logfc_thresh) {
      i <- i + 1
    } 
    else {
      gene_start <- num_conditions*i + 1
      gene_end <- num_conditions*(i + 1)
      de_cat_gene_df <- de_cat_gene_df[-c(gene_start:gene_end),]
    }
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
  }
  de_cat_gene_df$gene_timepoints <- rep(experimental_conditions, curr_num_genes)
 ggplot(de_cat_gene_df, aes(x=gene_timepoints, y=logFC, color=gene_name)) +
  scale_colour_brewer(palette="Paired")+
  geom_line(aes(color = gene_name), size = 1.25) +
  geom_point(aes(shape=below_fdr_threshold, color = gene_name), size = 5) +
  scale_shape_manual(values = c(1,17)) +
  scale_x_discrete(limits=experimental_conditions) +
  theme_bw() +
  ggtitle(plot_title) +
  theme(plot.title = element_text(hjust = 0.5, size = 36),
       legend.title = element_text(size = 23),
       legend.text = element_text(size =24),
       legend.position = "bottom",
       legend.box = "vertical",
       axis.text = element_text(size = 24),
       axis.title = element_text(size = 24)) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  coord_cartesian(ylim = plot_scale, xlim = c(0.25, 10.5)) +
  labs(shape = "Below FDR threshold", color = "Gene symbol") +
  xlab("Time post-infection [days]") + ylab("log2(Fold change)") +
  geom_hline(yintercept = 0)
}

gg1 <- plot_de_cat_genes2(cat_gene_list_p_u2, c(1,2,3,4,10), "IFN genes upregulated, primary infection", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-4,6))
gg2 <- plot_de_cat_genes2(cat_gene_list_s_u2, c(1,2,3,4,10), "IFN genes upregulated, secondary infection", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                  plot_scale = c(-4,6))
gg3 <- plot_de_cat_genes2(cat_gene_list_p_u1, c(1,2,3,4,10), "IFN genes upregulated, primary infection", fdr_thresh = fdr_threshold, 
                          logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                          plot_scale = c(-4,6))
gg4 <- plot_de_cat_genes2(cat_gene_list_s_u1, c(1,2,3,4,10), "IFN genes upregulated, secondary infection", fdr_thresh = fdr_threshold, 
                          logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                          plot_scale = c(-4,6))
gg5 <- plot_de_cat_genes2(cat_gene_list_p_d, c(1,2,3,4,10), "IFN genes downregulated, primary infection", fdr_thresh = fdr_threshold, 
                          logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                          plot_scale = c(-6,7))
gg6 <- plot_de_cat_genes2(cat_gene_list_s_d, c(1,2,3,4,10), "IFN genes downregulated, secondary infection", fdr_thresh = fdr_threshold, 
                          logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                          plot_scale = c(-6,7))
gg7 <- plot_de_cat_genes2(cat_gene_list_p_cc, c(1,2,3,4,10), "Cytokines & Chemokines, primary infection", fdr_thresh = fdr_threshold, 
                          logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                          plot_scale = c(-2,5))
gg8 <- plot_de_cat_genes2(cat_gene_list_s_cc, c(1,2,3,4,10), "Cytokines & Chemokines, secondary infection", fdr_thresh = fdr_threshold, 
                          logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0,
                          plot_scale = c(-2,5))

cat_plot_path <- "plots/old_and_new/mixed_IFN_up_expr_12.png"
png(cat_plot_path, width = 2000, height = 2000)
par(mar = c(5, 5, 5, 10))
ggarrange(ggarrange(gg2, gg1, labels = c("A", "B"),
                    font.label = list(size = 40),
                    ncol = 2, nrow = 1,
                    common.legend = TRUE,
                    legend="bottom"),
          ggarrange(gg4, gg3, labels = c("C", "D"),
                    font.label = list(size = 40),
                    ncol = 2, nrow = 1,
                    common.legend = TRUE,
                    legend="bottom"),
          ncol = 1, nrow = 2)
dev.off()
cat_plot_path <- "plots/old_and_new/mixed_IFN_down_expr.png"
png(cat_plot_path, width = 2000, height = 1000)
par(mar = c(5, 5, 5, 10))
ggarrange(gg6, gg5, 
          labels = c("A", "B"),
          font.label = list(size = 40),
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend="bottom")
dev.off()
cat_plot_path <- "plots/old_and_new/mixed_cyto_chemo_expr.png"
png(cat_plot_path, width = 2000, height = 1000)
par(mar = c(5, 5, 5, 10))
ggarrange(gg8, gg7,
          labels = c("A", "B"),
          font.label = list(size = 40),
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend="bottom")
dev.off()


# Plot volcanos side by side for primary & secondary:
volcano_path <- "plots/old_and_new/all_volcanos.png"
png(volcano_path, height = 1200, width = 2400)
ggarrange(ggarrange(plot_list_primary[[1]], plot_list_primary[[2]],plot_list_primary[[3]],
                    plot_list_primary[[4]],plot_list_primary[[5]],ncol = 2, nrow = 3),
          ggarrange(plot_list_secondary[[1]],plot_list_secondary[[2]],plot_list_secondary[[3]],
                    plot_list_secondary[[4]],plot_list_secondary[[5]],ncol = 2, nrow = 3),
          labels = c("A", "B"),
          font.label = list(size = 46),
          ncol = 2, nrow = 1)
dev.off()

# Plot volcanos day by day, side by side for primary & secondary:
volcano_path <- "plots/old_and_new/all_volcanos2.png"
png(volcano_path, height = 1900, width = 1400)
ggarrange(ggarrange(plot_list_primary[[1]], plot_list_primary[[2]],plot_list_primary[[3]],
                    plot_list_primary[[4]],plot_list_primary[[5]],ncol = 1, nrow = 5),
          ggarrange(plot_list_secondary[[1]],plot_list_secondary[[2]],plot_list_secondary[[3]],
                    plot_list_secondary[[4]],plot_list_secondary[[5]],ncol = 1, nrow = 5),
          labels = c("A", "B"),
          font.label = list(size = 46),
          ncol = 2, nrow = 1)
dev.off()
# Plot PCA side by side for primary+ secondary vs seconday infection:
pca_path <- "plots/old_and_new/mixed_pca.png"
pca_plot_p <- ggbiplot(pca_p, var.axes = FALSE, choices = 1:2, alpha = 1) +
  geom_text(aes(label = label_p), hjust = 0.5, vjust = -0.5, size = 4.5) +
  ggtitle("PCA of normalized chicken read counts per sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
pca_plot_s <- ggbiplot(pca_s, var.axes = FALSE, choices = 1:2, alpha = 1) +
  geom_text(aes(label = label_s), hjust = 0.5, vjust = -0.5, size = 4.5) +
  ggtitle("PCA of normalized chicken read counts per sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
png(pca_path, height = 1200, width = 2400)
ggarrange(pca_plot_p, pca_plot_s,
          labels = c("A", "B"),
          font.label = list(size = 46),
          ncol = 2, nrow = 1)
dev.off()

# WGCNA:
options(stringsAsFactors = FALSE)
chicken_expr <- t(cpm(chicken_dgelist_filt_norm))

# Pick out only genes from primary infection:
primary <- metadata[metadata$Infection_status=="primary",]
m <- match(primary$File_name, rownames(chicken_expr))
primary_chicken_expr <- chicken_expr[m,]

gsg = goodSamplesGenes(primary_chicken_expr, verbose = 3)
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
sampleTree = hclust(dist(primary_chicken_expr), method = "average")
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
chicken_expr_cut = primary_chicken_expr[keepSamples, ]
nGenes = ncol(chicken_expr_cut)
nSamples = nrow(chicken_expr_cut)

datExpr <- chicken_expr_cut
simMatrix <- cor(datExpr)
names(datExpr) = colnames(chicken_expr_cut);
rownames(datExpr) = rownames(chicken_expr_cut)
rm(chicken_expr_cut)
rm(chicken_cpm)
rm(simMatrix)


# Find soft treshold:
powers = c(c(1:10), seq(from = 11, to=20, by=1))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


# Plot soft treshold scale independence:
plot_path <- "plots/old_and_new/old_WGCNA_soft_tresh.png"
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
plot_path <- "plots/old_and_new/old_WGCNA_mean_con.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")
# Zoomed in
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"), ylim=range(0,500))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")
dev.off()


cor <- WGCNA::cor
# Run blockwise modules:
net9 <- blockwiseModules(datExpr, power = 9,
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE,
                          #saveTOMs = TRUE,
                          #saveTOMFileBase = "fusionTOM9", 
                          verbose = 3,
                          maxBlockSize = 18000)



table(net9$colors)
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net9$colors)


# Plot the dendrogram and the module colors underneath:
plot_path <- "plots/old_and_new/old_WGCNA_modules_block_net.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
plotDendroAndColors(net9$dendrograms[[1]], mergedColors[net9$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net9$colors
moduleColors = labels2colors(net9$colors)
dev.off()


MEs = net9$MEs
geneTree = net9$dendrograms[[1]]


# Set up trait file:
annotation <- metadata
rownames(annotation) <- annotation$File_name
dim(annotation)
names(annotation)
allTraits = annotation[, -c(3,4,5)];
dim(allTraits)
names(allTraits)
traitRows = match(rownames(datExpr), allTraits$File_name);
datTraits = allTraits[traitRows, -1, drop = FALSE];
rownames(datTraits) = allTraits[traitRows, 1];
traitRows = match(rownames(datExpr), allTraits$File_name);
datTraits = allTraits[traitRows, -1, drop = FALSE];
rownames(datTraits) = allTraits[traitRows, 1];


# Make binary matrix of timepoints:
datTraits[c(2,3,4,5,6)] <- datTraits[,1]
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


# Recalculate MEs with color labels:
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p", drop=FALSE)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)


# Will display correlations and their p-values:
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");


# Display the correlation values within a heatmap plot:
plot_path <- "plots/old_and_new/old_WGCNA_modules_time_p.png"
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
               main = paste("Module-trait relationships, primary infection"))
dev.off()

# Set up annotation file for only timepoint correlation:
annotation <- metadata
rownames(annotation) <- annotation$File_name
dim(annotation)
names(annotation)
allTraits = annotation[, -c(2, 4, 5)];
dim(allTraits)
colnames(allTraits) <- c( "File_name", "Timepoint")
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


# Create the starting data frame:
geneInfo0 <- data.frame(geneSymbol = colnames(datExpr),
                          EntrezID = chicken_dgelist_filt_norm$genes$entrez_gene_id[],
                          moduleColor = moduleColors,
                          geneTraitSignificance,
                          GSPvalue)


# Order modules by their significance for time_point:
modOrder <- order(-abs(cor(MEs, datTraits_time, use = "p")));
  

# Add module membership information in the chosen order:
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
  
write.csv(geneInfo, file = "tables/old_and_new/geneInfo_fusion.csv", row.names = FALSE)
  
geneInfo <- read.delim("tables/old_and_new/geneInfo_fusion.csv", sep = ",", stringsAsFactors = FALSE)
gene_info_mod <- geneInfo
chicken_product_df <- read.delim("chicken_gene_products.tsv", 
                                   header = FALSE, col.names = c("gene_name", "product"))
chicken_product_df <- chicken_product_df[!duplicated(chicken_product_df[,1]),]
  
egGENENAME <- toTable(org.Gg.egGENENAME)
gene_info_mod$Gene_name <- egGENENAME[match(gene_info_mod$EntrezID, egGENENAME$gene_id),]$gene_name
gene_info_mod$Product <- chicken_product_df[match(gene_info_mod$geneSymbol, chicken_product_df$gene_name),]$product
gene_info_mod$Differentially_expressed <- gene_info_mod$geneSymbol %in% logfc_filt_de_genes_chicken$gene_name
m <- match(gene_info_mod$geneSymbol, rownames(logfc_filt_de_genes_chicken))
gene_info_mod[,37:43] <- logfc_filt_de_genes_chicken[m,]
  
gene_info_mod <- gene_info_mod[,c(1,2,3,34:43,4:33)]
  
write.csv(gene_info_mod, file = "tables/old_and_new/geneInfo_fusion_annotated.csv", row.names = FALSE)
  
modules_genes <- read.delim("tables/old_and_new/geneInfo_fusion_annotated.csv", sep = ",", stringsAsFactors = FALSE)
DE_module_genes <- modules_genes[modules_genes$Differentially_expressed==TRUE ,]
write.csv(DE_module_genes, file = "tables/old_and_new/geneInfo_fusion_DE.csv", row.names = FALSE)
  
# GO and KEGG analyses of the modules
  
GOenr = GOenrichmentAnalysis(moduleColors, gene_info_mod$entrez_gene_id, organism = "chicken", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
write.table(tab, file = "tables/old_and_new/GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)
  
  
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
  # annotation should contain gene symbols and  either GO category or KEGG pathway membershp for both chicken and E. tenella genes  # cat_id is a string that tells what type of categories are being analysed
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
    go_file_name <- paste("tables/wgcna_module_go_kegg/old_and_new/", modNames[i], "_module_go_cats.csv", sep = "")
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
  kegg_file_name <- paste("tables/wgcna_module_go_kegg/old_and_new/", modNames[i], "_module_kegg_cats.csv", sep = "")
  write.csv(module_kegg_cats[[i]], file = kegg_file_name, row.names = FALSE)
}


  
  
modcol_old <- gene_types$module_old
modcol_new <- gene_types$module_col
counts1 <- table(modcol_old)
counts2 <- table(modcol_new)
bar_path <- "plots/old_and_new/IFN_module_bar.png"
png(bar_path, height = 800, width = 1700)
par(mfrow = c(1, 2))
par(mar = c(5, 5, 5, 5))
barplot(counts1, main="No. IFN-regulated genes in modules, primary",
        xlab="Module color", ylab="No. genes in module", ylim = c(0,20),
        col = c("blue", "brown","grey", "magenta", "purple", "red", "tan", "turquoise", "yellow"),
        cex.main=1.5, cex.lab=1.5, cex.axis=1)
barplot(counts2, main="No. IFN-regulated genes in modules, secondary",
        xlab="Module color", ylab="No. genes in module", ylim = c(0,20),
        col = c("blue", "brown", "cyan", "green", "grey", "lightcyan","red", "tan", "turquoise"),
        cex.main=1.5, cex.lab=1.5, cex.axis=1)
dev.off()
  