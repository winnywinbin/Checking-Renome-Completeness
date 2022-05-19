#install.packages("ggplot2", lib = "/usr/lib/R/library", type = 'source', dependencies=TRUE, repos = "http://cran.us.r-project.org")
library(ggplot2)


# Read file
#setwd("C:/Users/judyh/Documents/BPEXA")
mytable <- read.delim(snakemake@input[["filtered_archaea_data"]], header = TRUE, sep = "\t")
#mytable <- read.delim(snakemake@input[["domains_combined"]], header = TRUE, sep = "\t")
#mytable

# Alter table with 1 for presence anti-codon and 0 for absence
anticodons <- data.frame(unclass(table(mytable$GenomeID,mytable$Isotype)))
anticodons[anticodons > 1] <- 1
#anticodons

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

# Get sum of anti-codons for every column and create data frame
numbers <- c()
for (i in 1:ncol(anticodons)){
  numbers <- c(numbers, sum(anticodons[,i]))
}
df <- data.frame(
  anticodon=colnames(anticodons),
  n=numbers
)

# Create barplot and save as png
myplot <- ggplot(df, aes(x=anticodon, y=n)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90, hjust=1))
myplot
png(snakemake@output[["histogram"]])
print(myplot)
dev.off()