#debug(utils.get.unpackPkgZip)
#install.packages("rJava", lib = "/usr/lib/R/library", type = 'source', dependencies=TRUE, repos = "http://cran.us.r-project.org")
library(rJava)

install.packages("xlsx", dependencies = TRUE, repos = "http://cran.us.r-project.org")
#install.packages("contrib.url", repos = "http://cran.us.r-project.org")
library(xlsx)


# Read file
#setwd("C:/Users/judyh/Documents/BPEXA")
new_data <- read.delim(snakemake@input[["data"]], header = FALSE, sep = "\t")


# Alter table with 1 for presence anti-codon and 0 for absence
anticodons <- data.frame(unclass(table(new_data$V2,new_data$V5)))
anticodons[anticodons > 1] <- 1


# Hamming distance between genomes
# Function for calculating hamming distance
hamming.distance <- function(x,y){
  z<-NULL
  if(is.vector(x) && is.vector(y)){
    z <- sum(x != y)
  }
  else{
    z <- matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(k in 1:(nrow(x)-1)){
      for(l in (k+1):nrow(x)){
        z[k,l] <- hamming.distance(x[k,], x[l,])
        z[l,k] <- z[k,l]
      }
    }
    dimnames(z) <- list(dimnames(x)[[1]], dimnames(x)[[1]])
  }
  z
}


# Create empty matrix
number <- nrow(anticodons)
hmatrix <- matrix(, nrow = number, ncol=number)


# Calculate hamming distance between genomes and add to matrix
for (i in 1:nrow(anticodons)){
  G <- c()
  for (j in 1:nrow(anticodons)){
    i_vec <- as.numeric(as.vector(anticodons[i,]))
    j_vec <- as.numeric(as.vector(anticodons[j,]))
    G <- c(G, hamming.distance(i_vec, j_vec))
  }
  hmatrix[,i] <- G
}


# Rename rows and columns to the genome IDs
rownames(hmatrix) <- (rownames(anticodons))
colnames(hmatrix) <- (rownames(anticodons))


# Write to Excel file
write.xlsx(hmatrix, snakemake@output[["ha_excel"]])

