#setwd("C:/Users/daanv/PycharmProjects/project-bpexa")
#install.packages("gplots",lib = "/usr/lib/R/library", type = 'source', dependencies=TRUE, repos = "http://cran.us.r-project.org")
#install.packages("heatmaply",lib = "/usr/lib/R/library", type = 'source', dependencies=TRUE, repos = "http://cran.us.r-project.org")
library(gplots)
library(heatmaply)

archaea <- read.delim(snakemake@input[["filtered_archaea_data"]])
#mytable <- read.delim(snakemake@input[["filtered_archaea_data"]], header = TRUE, sep = "\t")
anticodons <- data.frame(unclass(table(archaea$Genome,archaea$Isotype)))

#Using gplots/heatmap.2

matrix <- as.matrix(anticodons)

heatmap1 <- heatmap.2(matrix)
htmlwidgets::saveWidget(as_widget(heatmap1), "test.html")

#Using heatmaply

df <- normalize(anticodons)
df[df == "NaN"] <- 0

#hm <- heatmaply(anticodons)
hmNormalized <- heatmaply(df, file="heatmap.html")
#htmlwidgets::saveWidget(as_widget(nmds_plot_i), snakemake@output[["nmdsplot"]])
#htmlwidgets::saveWidget(as_widget(hm), snakemake@output[["HeatmapAnticodons"]])
#htmlwidgets::saveWidget(as_widget(hmNormalized), snakemake@output[["HeatmapAnticodonsNormalized"]])
