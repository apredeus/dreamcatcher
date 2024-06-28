#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  cat(sprintf("Usage: ./parse_read_network.R <gcf-species-genus table> <component cutoff>\n"))
  quit(status=1) 
}

library(igraph)
library(dplyr)

## new version expects a network file in which there are self-edges

gsa <- read.csv(args[1], sep='\t', quote='', header=F, row.names=1)
colnames(gsa) <- c('Species_taxid','Species','Genus_taxid','Genus')
cutoff <- as.numeric(args[2])

filt_all <- try(read.csv('filtered.summary.tsv', sep='\t', header=F, quote='', comment.char = ''))
if (is.null(nrow(filt_all))) { 
  cat(sprintf("ERROR: no strains in filtered.summary.tsv! Exiting..\n"))
  quit(status=1)
}
colnames(filt_all) <- c('RefSeq','Strain_taxid','Strain','rRNA_genes','rRNA_reads','prot_genes','prot_reads')

edges_all <- try(read.table('nodes_and_edges.tsv', sep='\t'))
edges_strong <- edges_all[edges_all$V6 >= cutoff,] ## key parameter for this; currently set at 0.6
edges_strong <- edges_strong[edges_strong$V1 != edges_strong$V2,]
edges_strong <- edges_strong[,c(1,2,6)]
net <- graph_from_data_frame(d=edges_strong, directed=F)
comp <- components(net)

## in this version, all of the filtered strains are here obligatory
## here we get the table of strains to "honest reads" 
all_strains <- edges_all[,c(1,3)]
colnames(all_strains) <- c('RefSeq','Reads')
all_strains <- distinct(all_strains)

## now make a big table of all per-strain metadata, including the network 
all_strains <- merge(filt_all, all_strains, by='RefSeq')
all_strains <- merge(all_strains, as.data.frame(comp$membership), by.x='RefSeq', by.y=0, all.x=T)
colnames(all_strains)[9] <- 'Component'
all_strains[is.na(all_strains)] <- 'N'
all_strains <- merge(all_strains, as.data.frame(comp$csize), by.x='Component', by.y=0, all.x=T)
all_strains[is.na(all_strains)] <- '1'
colnames(all_strains)[10] <- 'Component_size'

## split into 2 parts - 1st will be collapsed, 2nd used as is
net_strains  <- all_strains[all_strains$Component != 'N',]
uniq_strains <- all_strains[all_strains$Component == 'N',]

net_strains <- net_strains[order(net_strains$Component, -net_strains$Reads),]
net_strains <- net_strains[!duplicated(net_strains$Component),]
top_strains <- rbind(net_strains, uniq_strains)
top_strains <- merge(top_strains, gsa, by.x='RefSeq', by.y=0)
top_strains$Component_size <- as.numeric(top_strains$Component_size)
top_strains <- top_strains[with(top_strains, order(-Component_size,-rRNA_reads)), ]

rownames(top_strains) <- 1:nrow(top_strains)
top_strains <- top_strains[,c(1,11:14,5:9,2,10)]

write.table(top_strains,'top_strains.tsv', quote=F, sep='\t', row.names=F)

