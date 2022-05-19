# Read file
#setwd("C:/Users/judyh/Documents/BPEXA")
#mytable <- read.delim(snakemake@input[["filtered_archaea_data"]], header = TRUE, sep = "\t")
mytable <- read.delim("Filtered_archaea_GtRNAdb.txt", header = TRUE, sep = "\t")

# Alter table with 1 for presence anti-codon and 0 for absence
anticodons <- data.frame(unclass(table(mytable$GenomeID,mytable$Isotype)))
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

# Get sum of anti-codons for every column and create data frame
numbers <- c()
for (i in 1:ncol(anticodons)){
  numbers <- c(numbers, sum(anticodons[,i]))
}
df <- data.frame(
  anticodon=colnames(anticodons),
  n=numbers
)

# Create list of expected anticodons and list of not expected anticodons
hpercentage <- as.integer(nrow(anticodons)/100*95)
lpercentage <- as.integer(nrow(anticodons)/100*5)
not.expected <- c()
expected <- c()
for (i in 1:nrow(df)){
  number <- df$n[i]
  if (number >= hpercentage){
    expected <- c(expected,df$anticodon[i])
  }
  else if (number <= lpercentage){
    not.expected <- c(not.expected,df$anticodon[i])
  }
}

# Check completeness
# Can only use files with same domain as data in dataframe 'mytable'
#genomes <- read.delim(snakemake@input[["filtered_data"]],header=FALSE,sep="\t")
genomes <- read.delim("FilteredData.txt",header=FALSE,sep="\t")
ac <- data.frame(unclass(table(genomes$V2,genomes$V6)))
ac[ac > 1] <- 1

# Loop through every row (genome) in dataframe to calculate 
# completeness and write to file
for (i in 1:nrow(ac)){
  genomeid <- row.names(ac[i,])
  list <- colnames(ac[i,])[ac[i,] == 1]
  count=0
  l.not.expected <- c()
  
  # Check if anticodons are expected or not
  for (a in list){
    if (a %in% expected){
      count=count+1
    }
    else if (a %in% not.expected){
      l.not.expected <- c(l.not.expected,a)
    }
  }
  
  # Check if anticodons are missing
  left <- c()
  for (a in expected){
    if (a %in% list){
      next
    } else {
      left <- c(left,a)
    }
  }
  # Calculate completeness
  comp <- round(count/length(expected)*100,digits=1)
  
  #Write to file
  myfile <- snakemake@output[["completeness"]]
  c <- paste0("Completeness: ",comp,"%")
  cat(genomeid,c,"Expected but absent anticodons: ",left,"Not expected but present anticodons: ",
      l.not.expected," ",file=myfile,sep="\n",append=TRUE)
}

