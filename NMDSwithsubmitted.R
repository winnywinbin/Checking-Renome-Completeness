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

mytable <- read.delim(snakemake@input[["domains_combined"]], header = TRUE, sep = "\t")
genomes <- read.delim(snakemake@input[["data"]], header = FALSE, sep = "\t")
mytable$ID <- paste(mytable$GenomeID, " (", mytable$Domain ,")")
genomes$ID <- paste(genomes$V2, " (", "Submitted",")")

# Alter table with 1 for presence anti-codon and 0 for absence
anticodons <- data.frame(unclass(table(mytable$ID,mytable$Isotype)))
anticodons.genomes <- data.frame(unclass(table(genomes$ID,genomes$V6)))
anticodons[anticodons > 1] <- 1
anticodons.genomes[anticodons.genomes > 1] <- 1

# Add absent anticodons to matrices and combine
allanticodons <- c('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
                   'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
                   'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
                   'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
                   'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
                   'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
                   'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
                   'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT')

add_columns = function(matrix){
  zeros <- rep(c(0), times=nrow(matrix))
  for (anticodon in allanticodons){
    if (anticodon %in% colnames(matrix)){
      next
    } else {
      matrix <- cbind(matrix, zeros)
      names(matrix)[ncol(matrix)] <- anticodon
    }
  }
  return(matrix)
}

anticodons <- add_columns(anticodons)
anticodons.genomes <- add_columns(anticodons.genomes)
anticodons <- rbind(anticodons, anticodons.genomes)

# NMDS plot
nmds_results <- metaMDS(comm=anticodons, distance="bray", k=8, try=10)

data_scores <- as.data.frame(scores(nmds_results))
data_scores$Domain <- sapply(strsplit(rownames(data_scores), " "), "[", 4)
data_scores$GenomeID <- sapply(strsplit(rownames(data_scores), " "), "[", 1)

nmds_plot <- ggplot(data=data_scores, aes(x=NMDS1, y=NMDS2, label=GenomeID, color=Domain)) +
  geom_point(aes(shape=Domain)) + 
  scale_shape_manual(values=c(1, 4, 17, 20)) +
  ggtitle("NMDS plot of genomes")
nmds_plot_i <- ggplotly(nmds_plot)
nmds_plot_i


htmlwidgets::saveWidget(as_widget(nmds_plot_i), "nmdsplotwithsubmitted.html")






