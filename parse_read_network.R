library(igraph)
library(dplyr)

tt <- read.table('nodes_and_edges.tsv', sep='\t')
ttf <- tt[tt$V6 >= 0.5,]
ttw <- ttf[,c(1,2,6)]
net <- graph_from_data_frame(d=ttw, directed=F)
comp <- components(net)

strains1 <- ttf[,c(1,3)]
strains2 <- ttf[,c(2,4)]
colnames(strains1) <- c('Strain','Reads')
colnames(strains2) <- c('Strain','Reads')
strains <- distinct(rbind(strains1,strains2))
dim(strains)
rownames(strains) <- strains$Strain
strains <- merge(strains,as.data.frame(comp$membership),by=0)
strains$Row.names <- NULL
colnames(strains)[3] <- 'Component'
strains <- merge(strains,as.data.frame(comp$csize),by.x='Component',by.y=0)


##sort and find the biggest guy
ordered_strains <- strains[order(strains$Component, -strains$Reads),]
top_strains <- ordered_strains[!duplicated(ordered_strains$Component),]
colnames(top_strains)[4] <- 'Component_size'
write.table(top_strains,"top_strains.tsv", quote=F, sep='\t', row.names=F)
