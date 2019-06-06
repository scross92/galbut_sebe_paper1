#Correlation analyses on 2017/2018 Drosophila
#This data does not include individual Brazil flies


library(ggplot2)
library(dplyr)
library(googlesheets)
library(psych)

# open google spreadsheet
flies_gs <- gs_title("drosophila_colony_taxa_matrix")
gs_ws_ls(partiti_gs)
fly_taxa_matrix <- gs_read(ss=flies_gs)

##The following removes host reads, and then collapses reads based on genus name
##Probably should have done this in the terminal rather than here...just to have it cleaner

#remove columns that are host reads
fly_taxa_matrix_sub <- select(fly_taxa_matrix, -contains("Drosophila "))
#View(fly_taxa_matrix_sub)
fly_taxa_matrix_sub[is.na(fly_taxa_matrix_sub)] <- 0

#remove all remaining Wolbachia reads
#Most wolbachia get removed with host name pattern
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Wolbachia"))

#Add all Wolbachia columns
#Need to do this separately because most of the Wolbachia have Drosophila in their name, so they get removed when removing host reads using
#the pattern "Drosophila"
wolbachia_matrix_sub <- select(fly_taxa_matrix, contains("Wolbachia"))
col_number <- ncol(wolbachia_matrix_sub)
wolbachia_matrix_sub[is.na(wolbachia_matrix_sub)] <- 0
#View(wolbachia_matrix_sub)

patterns <- unique(substr(names(wolbachia_matrix_sub), 1 , 6))  # store patterns in a vector
wolbachia_matrix_sub_sum <- sapply(patterns, function(xx) rowSums(wolbachia_matrix_sub[,grep(xx, names(wolbachia_matrix_sub)), drop=FALSE])) #loop through
###colnames(new) <- paste0(colnames(new), "tot")  # rename by adding tot to the end of each column
#View(wolbachia_matrix_sub_sum)

#Add all Acetobacter
acetobacter_sums <- rowSums(select(fly_taxa_matrix_sub, contains("Acetobacter")))

#Add all Lactobacillus
lactobacillus_sums <- rowSums(select(fly_taxa_matrix_sub, contains("Lactobacillus")))

#Add all Saccharomyces
saccharomyces_sums <- rowSums(select(fly_taxa_matrix_sub, contains("Saccharomyces")))

#Add all Gluconobacter
gluconobacter_sums <- rowSums(select(fly_taxa_matrix_sub, contains("Gluconobacter")))

#Add all Hanseniaspora
hanseniaspora_sums <- rowSums(select(fly_taxa_matrix_sub, contains("Hanseniaspora")))

#Add all Lactococcus
lactococcus_sums <- rowSums(select(fly_taxa_matrix_sub, contains("Lactococcus")))

#Remove all multiples
#this does it one at a time...not sure if there is a way to make this faster...
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Acetobacter"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Lactobacillus"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Saccharomyces"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Gluconobacter"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Hanseniaspora"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Lactococcus"))
#View(fly_taxa_matrix_sub)

#Remove "uncultured" species
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("uncultured "))
#View(fly_taxa_matrix_sub)
#Remove Unknown species and generic bacterium
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Unknown"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("bacterium"))
#View(fly_taxa_matrix_sub)

#Add back in cols that were removed
fly_taxa_matrix_sub <- cbind(fly_taxa_matrix_sub, acetobacter_sums, lactobacillus_sums, saccharomyces_sums, gluconobacter_sums, hanseniaspora_sums, lactococcus_sums, wolbachia_matrix_sub_sum)

#Remove nonsense human herpesvirus, and other garbage reads (carp, flowers, etc)
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Human"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Cyprinus carpio"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Glyptapanteles flavicoxis"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Lasthenia californica"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Avena fatua"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Phleum pratense"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Brachypodium distachyon"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Panicum hallii"))
fly_taxa_matrix_sub <- select(fly_taxa_matrix_sub, -contains("Penaeus vannamei"))



#remove names in the first column
fly_taxa_matrix_sub_nameless <- fly_taxa_matrix_sub[,-1]

#normalize the data by unique reads per million
fuh_reads <- read.csv(file.choose(), row.names = 1)

#remove rows that do not contain any values (they didn't have any nonhost reads)
fly_taxa_matrix_sample_subset <- fly_taxa_matrix_sub_nameless[apply(fly_taxa_matrix_sub_nameless[,-1], 1, function(x) !all(x==0)),]
#View(fly_taxa_matrix_sample_subset)

fly_taxa_matrix_norm <- (fly_taxa_matrix_sample_subset / fuh_reads[col(fly_taxa_matrix_sample_subset)]) * 1e6
View(fly_taxa_matrix_norm)

#remove virus reads
fly_taxa_matrix_norm_no_viruses <- select(fly_taxa_matrix_sub, -contains("virus"))
fly_taxa_matrix_norm_no_viruses <- select(fly_taxa_matrix_norm_no_viruses, -contains("Partiti"))

#remove fungal reads (bacteria only)
fly_taxa_matrix_bacteria <- select(fly_taxa_matrix_norm_no_viruses, -contains("saccharomyces"))
fly_taxa_matrix_bacteria <- select(fly_taxa_matrix_bacteria, -contains("hanseniaspora"))
fly_taxa_matrix_bacteria <- fly_taxa_matrix_bacteria[,-1]

#calculate total number of bacterial mapping reads per sample
bacteria_reads <- rowSums(fly_taxa_matrix_bacteria)
#calculate avg # of bacteria reads per sample
bacteria_reads_avg <- (sum(bacteria_reads))/(NROW(bacteria_reads))
bacteria_reads_avg
#calculate median # of bacteria reads
median(bacteria_reads)

#virus reads only matrix
fly_taxa_matrix_viruses <- select(fly_taxa_matrix_sub, contains("virus"))
partiti2 <- select(fly_taxa_matrix_sub, contains("Partitiviridae"))
fly_taxa_matrix_viruses <- cbind(fly_taxa_matrix_viruses, partiti2)

#calculate total number of virus mapping reads per sample
virus_reads <- rowSums(fly_taxa_matrix_viruses)
#calculate avg # of virus reads per sample
virus_reads_avg <- (sum(virus_reads))/(NROW(virus_reads))
virus_reads_avg

#bacteria vs virus reads scatter plot with basic plot function
bacteria_virus_reads <- data.frame(bacteria_reads, virus_reads)
# plot(bacteria_virus_reads$bacteria_reads, bacteria_virus_reads$virus_reads, 
#      xlab = "Bacterial Reads per Million Reads", 
#      ylab = "Viral Reads Per Million Reads", 
#      main = "Bacteria vs Virus Mapping Reads in Individual Flies",
#      yaxt = "n",
#      xaxt = "n")
# 
# #Fix the axis to have better labels
# # Define the position of tick marks
# y1 <- c(0,2e4,4e4,6e4,8e4, 1e5)
# 
# # Define the labels of tick marks
# y2 <- c("0","20,000","40,000","60,000","80,000", "100,000")
# 
# # Add an axis to the plot--ylab 
# axis(side = 2, 
#      at = y1, 
#      labels = y2,
#      tck=-.02)
# 
# #Position and labels
# x1 <- c(0, 5000, 10000, 15000, 2e4)
# x2 <- c("0", "5,000", "10,000", "15,000", "20,000")
# #Add an axis to plot--xlab
# axis(side = 1, 
#      at = x1, 
#      labels = x2,
#      tck=-.02)
# 
# #Add fitted line (maybe not the best?)
# reads_lm <- lm(bacteria_virus_reads$virus_reads~bacteria_virus_reads$bacteria_reads)
# abline(reads_lm)


#making a scatterplot with ggplot2
#ggplot2 has more versatility
bacteria_virus_reads
library(grid)

p <- ggplot(bacteria_virus_reads, aes(x = bacteria_reads, y = virus_reads)) + geom_point()
p + geom_abline(intercept = 0) + #add y=x line
  theme_bw() + #remove background color
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(color = "black")) + #remove grid lines and border lines
  labs(x = "Bacteria Reads per Million Reads", y = "Virus Reads per Million Reads") + #axis labels
  theme(axis.title = element_text(family = "Helvetica", size = 15)) +
  ggtitle("Bacteria vs Virus Mapping Reads in Individual Drosophila") + #title
  theme(plot.title = element_text(family = "Helvetica", size = 20, hjust = 0.5)) + #Center the title
  theme(axis.text = element_text(size = 12)) + #change axis text size
  theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5),"cm")) + #add margins around plot
  theme(
    axis.title.x = element_text(margin = unit(c(1, 0, 0.5, 0), "cm")), # add space between x label and plot
    axis.title.y = element_text(margin = unit(c(0, 1, 0, 0.5), "cm")), #add space between y lab and plot
    plot.title = element_text(margin = unit(c(0.5, 0, 1.5, 0), "cm"))) #add space between title and plot
  

#Corr.test from psych library returns pearson coefficients and p-vals

fly_taxa_matrix_norm[fly_taxa_matrix_norm == 0] <- NA
fly.corr.test <- corr.test(fly_taxa_matrix_norm, use = "pairwise", method = "pearson", adjust = "bonferroni")
write.csv(fly.corr.test$n,file="ntmp.csv")  ## sample sizes
write.csv(fly.corr.test$t,file="ttmp.csv")  ## t statistics
write.csv(fly.corr.test$p,file="ptmp.csv")  ## p-values
write.csv(fly.corr.test$r, file = "rtmp.csv") ## coefficients
