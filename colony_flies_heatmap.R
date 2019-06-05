#This script makes a heatmap of all the viruses in each individual 2017/2018 fly
#This does not include brazil fly data

library(RColorBrewer)
library(gplots)

# export normalized data
hm <- read.delim(file.choose(), sep="\t",row.names=1)
hm
# convert to matrix (not dataframe)
hm_mat <- data.matrix(hm)

# log transform values
hm_log <- log(hm_mat, base=10)
# zero out infinite values
hm_log [which(!is.finite(hm_log))] <- 0

# creates a blue color pallette
my_palette <- brewer.pal(7, "Blues")
# replace the first light blue color in this palette w/ actual white
my_palette <- replace(my_palette, 1, "#FFFFFF")

heatmap.2(hm_log,
          main = "Virome of individual wild Drosophila melanogaster", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color 
          legend,
          trace="none",         # turns off trace lines inside the heat map
          margins =c(8,25),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram='none',    # No dendogram
          Colv=FALSE,
          Rowv=FALSE,           # Don't cluster/reorder rows
          sepwidth=c(0.02,0.1),
          sepcolor="grey90",
          colsep=1:ncol(hm_log),
          rowsep=1:nrow(hm_log),
          cexCol = 0.9)       


