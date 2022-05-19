#install.packages("gganimate", lib = "/usr/lib/R/library", type = 'source', dependencies=TRUE, repos = "http://cran.us.r-project.org")
#install.packages("dplyr", lib = "/usr/lib/R/library", type = 'source',  dependencies=TRUE, repos = "http://cran.us.r-project.org")
#install.packages("ggplot2", lib = "/usr/lib/R/library", type = 'source',  dependencies=TRUE, repos = "http://cran.us.r-project.org")
library(dplyr)
library(ggplot2)
library(gganimate)

#setwd("~/Desktop/project-bpexa")

new_dataset <- read.delim(snakemake@input[["filtered_data"]], header = F,  fill = TRUE)

# construct dataframe with counts
tRNAcount.tmp <-
  new_dataset %>% group_by(V1) %>% summarise(nr_tRNAs = length(V1))

anticodoncount.tmp <-
  new_dataset %>% group_by(V1) %>% summarise(no_anticodons = n_distinct(V6))

data_counts <- merge(tRNAcount.tmp, anticodoncount.tmp, by="V1")

#data_counts <- data_counts.tmp[-c(1),]

write.table(data_counts, file = snakemake@output[["data_counts"]])

#---- Archaea histograms  

# number of tRNAs
ggplot(data_counts, aes(nr_tRNAs)) +
  geom_histogram() +
  scale_x_continuous(n.breaks=2) +
  ggtitle("tRNA count")

# number of anticodons
ggplot(data_counts, aes(no_anticodons)) +
  geom_histogram() +
  scale_x_continuous(n.breaks=20) +
  ggtitle("Anticodon count")

#---- Archaea filtering
# 36 species with few tRNAs
#data_counts %>% filter(nr_tRNAs < 46)
# 53 species with missing anticodons
#data_counts %>% filter(no_anticodons < 44)

new_dataset <- transform(new_dataset, count = ave(V4, V1, FUN = length))

data_anticodons <-
  new_dataset %>% group_by(V1, V6, count) %>% tally()

barplot <-
  ggplot(data_anticodons, aes(x=V6, y=n, fill=V6)) +
  geom_bar(stat = "summary", fun="mean") +
  theme(axis.text.x = element_text(angle=90),
        legend.position="none") + 
  labs(x="Anticodons", y="Count") +
  #transition_time(count) +
  labs(title = "tRNAs: {frame_time}")

#animate(barplot, duration=30, nframes=30, renderer=gifski_renderer("data_anticodon_plot.gif"))
#animate(barplot)

#---- Archaea fisher exact tests

# print only the relevant columns
data_twocol <- new_dataset[,c("V1","V6")]
# tabulate and save as dataframe
# 53 anticodons, 217 species
data_table <- data.frame(unclass(table(data_twocol)))

# change to presence/absence
data_table[data_table > 1] <- 1
#data_table <- rbind(data_table, "000" = c(0))
# fisher tests give error because most 2x2 tables are incomplete
# combn(archaea_table,2, function(x) fisher.test(table(x)), simplify=F)

# workaround of the above problem
# first calculate 2x2 table for each comination of anticodons
data_2x2 <- combn(data_table,2, function(x) table(x), simplify=F)

# then conduct fisher on each table if possible
sink(file = snakemake@output[["fisher_output_R"]])
for (i in data_2x2){
  tryCatch(
    expr={
      print(fisher.test(i))
      print(i)
    }, 
    error=function(e){print("NULL")})
}
sink()
closeAllConnections()