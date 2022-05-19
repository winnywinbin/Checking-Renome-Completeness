#install.packages("xlsx")
#install.packages("contrib.url")
library(xlsx)


# Read file
#setwd("~/Desktop/project-bpexa")
new_data <- read.delim(snakemake@input[["data"]], header = FALSE, sep = "\t")


# Alter table with 1 for presence anti-codon and 0 for absence
anticodons <- data.frame(unclass(table(new_data$V2,new_data$V6)))
anticodons[anticodons > 1] <- 1


# Cosine similarity between genomes
# Function for calculating cosine similarity
cos_sim = function(x, y){
  res = x %*% y / (sqrt(x %*% x) * sqrt(y %*% y))
  res = as.numeric(res)
  return(res)
}


# Create empty matrix
number <- nrow(anticodons)
gmatrix <- matrix(, nrow = number, ncol=number)


# Calculate cosine similarity between genomes and add to matrix
for (i in 1:nrow(anticodons)){
  G <- c()
  for (j in 1:nrow(anticodons)){
    i_vec <- as.numeric(as.vector(anticodons[i,]))
    j_vec <- as.numeric(as.vector(anticodons[j,]))
    G <- c(G, cos_sim(i_vec, j_vec))
  }
  gmatrix[,i] <- G
}


# Rename rows and columns to the genome IDs
rownames(gmatrix) <- (rownames(anticodons))
colnames(gmatrix) <- (rownames(anticodons))


# Write to Excel file
write.xlsx(gmatrix, snakemake@output[["co_excel"]])

