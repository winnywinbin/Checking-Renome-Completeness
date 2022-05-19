install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)

library(phytools)


#newick code inlezen
tree<-read.tree("domain.ph")

#data inlezen die tegen de boom moet worden geplot, moet matrix zijn anders kan de module er niet mee werken
#data met 1 datapunt
#data2<-as.matrix(read.csv("data.csv", row.names=1))[,1]
 
#data met twee datapunten
data<-read.delim("tdata_clean.csv", row.names=1, sep = ",")


#collapse tree
#collapsed <- di2multi(tree, tol = 0.5)
#dotTree(collapsed, data,length=10,ftype="i")
#pl<-plot(collapsed,type = "phylogram")

#boom wegschrijven als .png file 
png(filename = "tree_met_data.png", res = 100,
    width = 12000, height = 10000)
pl<-dotTree(tree,data,length=10,ftype="i")
dev.off()
