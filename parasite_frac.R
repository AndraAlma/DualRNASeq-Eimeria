
path_to_221221 <- "~/exjobb/metadata_table_short_221221.csv"
path_to_221111 <- "~/exjobb/metadata_table_short_221111.csv" 
path_to_221228 <- "~/exjobb/metadata_table_short_221228.csv" 
path_to_230119 <- "~/exjobb/metadata_table_short_230119.csv" 
path_to_pilot <- "~/exjobb/metadata_table_short_pilot.csv" 

#read in 
stats_221221 <- read.delim(path_to_221221, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
stats_221111 <- read.delim(path_to_221111, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
stats_221228 <- read.delim(path_to_221228, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
stats_230119 <- read.delim(path_to_230119, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
stats_pilot <- read.delim(path_to_pilot, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
stats_221221 <- merge(stats_221221,chicken_metadata, by.x = "File_name", by.y = "File_name")
stats_221111 <- merge(stats_221111,chicken_metadata, by.x = "File_name", by.y = "File_name")
stats_221228 <- merge(stats_221228,chicken_metadata, by.x = "File_name", by.y = "File_name")
stats_230119 <- merge(stats_230119,chicken_metadata, by.x = "File_name", by.y = "File_name")
stats_pilot <- merge(stats_pilot,chicken_metadata, by.x = "File_name", by.y = "File_name")

tot_stats <- rbind(stats_221221, stats_221111, stats_221228, stats_230119, stats_pilot)

plot(tot_stats$Timepoint, tot_stats$Percentage_of_Eimeria_reads, xaxt = "n", col= tot_stats$`Sample date`)
axis(1, at=0:4, labels=0,1,2,3,4,10)
df2 <- data.frame(x = tot_stats$Timepoint , y = tot_stats$Percentage_of_Eimeria_reads, z=tot_stats$`Sample date`)
df2$z <- as.factor(df2$z)



plot_path <- "~/Exjobb/plots/Eim_frac_bw.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
ggplot(df2, aes(x, y)) +
  geom_point(size = 3) +
  scale_x_break(c(4, 9.5))+
  xlim(0, 10.5)+
  scale_x_continuous(breaks = c(0,1,2,3,4,10),limits = c(0, 10.25), labels=c("Day -3","Day 1","Day 2","Day 3","Day 4","Day 10"))+
  xlab("Timepoint (Day)") + 
  ylab("% Eimeria reads")
dev.off()


plot_path <- "~/Exjobb/plots/Eim_frac_col.png"
png(plot_path, height = 600, width = 800, pointsize = 16)
ggplot(df2, aes(x, y, color=z)) +
  geom_point(size = 3) +
  scale_x_break(c(4, 9.5))+
  xlim(0, 10.5)+
  scale_x_continuous(breaks = c(0,1,2,3,4,10),limits = c(0, 10.25), labels=c("-3","1","2","3","4","10"))+
  xlab("Timepoint (Day)") + 
  ylab("% Eimeria reads")+
  scale_fill_brewer(palette="Dark2")
dev.off()
 
  
  

