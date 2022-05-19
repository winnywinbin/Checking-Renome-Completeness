#install.packages("ggplot2", lib = "/usr/lib/R/library", type = 'source', dependencies=TRUE, repos = "http://cran.us.r-project.org")
#install.packages("plotly", lib = "/usr/lib/R/library", type = 'source', dependencies=TRUE, repos = "http://cran.us.r-project.org")
#install.packages("vegan", lib = "/usr/lib/R/library", type = 'source', dependencies=TRUE, repos = "http://cran.us.r-project.org")
#install.packages("viridis", lib = "/usr/lib/R/library", type = 'source', dependencies=TRUE, repos = "http://cran.us.r-project.org")
#install.packages("ggplot2")
#install.packages("plotly")
#install.packages("vegan")
#install.packages("viridis")
library(ggplot2)
library(plotly)
library(vegan)
library(viridis)

# Read file
#setwd("C:/Users/judyh/Documents/BPEXA")
mytable <- read.delim(snakemake@input[["domains_combined"]], header = TRUE, sep = "\t")
#mytable <- read.delim(snakemake@input[["filtered_data"]], header = F,  fill = TRUE)
mytable$ID <- paste(mytable$GenomeID, " (", mytable$Domain ,")")

# Alter table with 1 for presence anti-codon and 0 for absence
anticodons <- data.frame(unclass(table(mytable$ID,mytable$Isotype)))
anticodons[anticodons > 1] <- 1

# Add absent anticodons to matrix
allanticodons <- c('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
                   'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
                   'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
                   'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
                   'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
                   'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
                   'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
                   'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT')
zeros <- rep(c(0), times=nrow(anticodons))
for (anticodon in allanticodons){
  if (anticodon %in% colnames(anticodons)){
    next
  } else {
    anticodons <- cbind(anticodons, zeros)
    names(anticodons)[ncol(anticodons)] <- anticodon
  }
}

# # PCA
# pca <- prcomp(anticodons) 
# pca.data <- data.frame(Sample=rownames(pca$x),X=pca$x[,2],Y=pca$x[,3])
# pca.data$Domain <- sapply(strsplit(pca.data$Sample, " "), "[", 4)
# 
# p <- ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, color=Domain)) +
#   geom_point(shape=1) + 
#   xlab("pca$x[,2]") + 
#   ylab("pca$x[,3]") + 
#   ggtitle("PCA of all genomes")
# pl <- ggplotly(p)
# htmlwidgets::saveWidget(as_widget(pl), "plot_23.html")
# 
# # Variance explained barplot
# p.variance.explained = pca$sdev^2 / sum(pca$sdev^2)
# # plot percentage of variance explained for each principal component    
# b <- barplot(100*p.variance.explained, las=2, names.arg=colnames(pca$rotation), ylab='% Variance Explained')


# NMDS plot
nmds_results <- metaMDS(comm=anticodons, distance="bray", try=10)

data_scores <- as.data.frame(scores(nmds_results))
data_scores$Domain <- sapply(strsplit(rownames(data_scores), " "), "[", 4)
data_scores$GenomeID <- sapply(strsplit(rownames(data_scores), " "), "[", 1)

nmds_plot <- ggplot(data=data_scores, aes(x=NMDS1, y=NMDS2, label=GenomeID, color=Domain)) +
  geom_point(aes(shape=Domain)) + 
  scale_shape_manual(values=c(1, 4, 17)) +
  ggtitle("NMDS plot of genomes")
nmds_plot_i <- ggplotly(nmds_plot)

htmlwidgets::saveWidget(as_widget(nmds_plot_i), snakemake@output[["nmdsplot"]])
#htmlwidgets::saveWidget(as_widget(nmds_plot_i), snakemake@output[["nmdsplot.html"]])
#write.table(data_counts, file = snakemake@output[["data_counts"]])

