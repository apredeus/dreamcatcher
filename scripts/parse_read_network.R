#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  cat(sprintf("Usage: ./parse_read_network.R <gcf_species_genus.tsv>\n"))
  quit(status=1)
}

library(igraph)
library(dplyr)

### update on the clustering vNNN:
##  - select all of the good strains - uniq_prot > 0 and prot_genes > 2
##  - only count the edges with weight of 0.5 or more
##  - cluster them with Leiden using absolute overlap as a weight
##  - get clusters that have edge_density of 0.9 or more, use as is
##  - split clusters that have edge_denisty of < 0.9, to subclusters == genus
##  - calculate score: uniq_rrna + 10*uniq_prot; take the strain with highest score as representative for each cluster

##  - group all of the remaining (non "uniq") strains into clusters == genus
##  - iterate through _all_ clusters seing if they can be combined and still retain edge_density of 0.9 or more
##  - stop when no pairs that can be combined in this way remain

### first let's make the big metadata table, 1 row = 1 filtered strain

gsa        <- try(read.csv(args[1], sep='\t', header=F, quote='', comment.char='', row.names=1))
filt_all   <- try(read.csv("filtered.summary.tsv", sep='\t', header=F, quote='', comment.char=''))
uniq_rrna  <- try(read.table("rRNA.uniq_counts.tsv", sep='\t'))         ## unique rRNA reads per strain
uniq_prot  <- try(read.table("protein.uniq_counts.tsv", sep='\t'))      ## unique non-rRNA reads per strain
edges_all  <- try(read.table("nodes_and_edges.tsv", sep='\t'))          ## all shared reads, including loops (strain1 == strain2)
edges_rrna <- try(read.table("rRNA.nodes_and_edges.tsv", sep='\t'))     ## all shared rRNA reads
edges_prot <- try(read.table("protein.nodes_and_edges.tsv", sep='\t'))  ## all shared non-rRNA reads

colnames(filt_all) <- c('RefSeq','Strain','Strain_taxid','rRNA_genes','rRNA_fcount','rRNA_raw','rRNA_mis','rRNA_hum',
                        'prot_genes','prot_fcount','prot_raw','prot_mis','prot_hum')

strain_ann <- filt_all[,c(1,4,5,9,10)]
colnames(gsa) <- c('Species_taxid','Species','Genus_taxid','Genus')

## check if any of the uniq files are empty, replace with NULL if so

if (any(grepl("error", uniq_rrna, ignore.case=T))) { 
  cat(sprintf("WARNING: No unique rRNA reads found for this experiment!\n"))
  uniq_rrna <- NULL
} else { 
  colnames(uniq_rrna) <- c("RefSeq", "uniq_rrna")
}
if (any(grepl("error", uniq_prot, ignore.case=T))) { 
  cat(sprintf("WARNING: No unique protein reads found for this experiment!\n"))
  uniq_prot <- NULL
} else { 
  colnames(uniq_prot) <- c("RefSeq", "uniq_prot")
}

## nodes_and_edges.tsv should always exist for a non-empty filtered strain list; other two, not so much

if (any(grepl("error", edges_rrna, ignore.case=T))) { 
  cat(sprintf("WARNING: No shared rRNA reads found for this experiment!\n"))
  shared_rrna <- NULL
} else { 
  shared_rrna <- edges_rrna[edges_rrna$V1 == edges_rrna$V2, c(1,3)]
  colnames(shared_rrna) <- c("RefSeq", "shared_rrna")
}
if (any(grepl("error", edges_prot, ignore.case=T))) { 
  cat(sprintf("WARNING: No shared protein reads found for this experiment!\n"))
  shared_prot <- NULL
} else { 
  shared_prot <- edges_prot[edges_prot$V1 == edges_prot$V2, c(1,3)]
  colnames(shared_prot) <- c("RefSeq", "shared_prot")
}

## merge all of those into the main ann
if (is.null(uniq_rrna)) { 
  strain_ann$uniq_rrna <- 0
} else {
  strain_ann <- merge(strain_ann, uniq_rrna, by="RefSeq", all=T)
}
if (is.null(uniq_prot)) { 
  strain_ann$uniq_prot <- 0
} else {
  strain_ann <- merge(strain_ann, uniq_prot, by="RefSeq", all=T)
}
if (is.null(shared_rrna)) { 
  strain_ann$shared_rrna <- 0
} else {
  strain_ann <- merge(strain_ann, shared_rrna, by="RefSeq", all=T)
}
if (is.null(shared_prot)) { 
  strain_ann$shared_prot <- 0
} else {
  strain_ann <- merge(strain_ann, shared_prot, by="RefSeq", all=T)
}

strain_ann <- merge(strain_ann, gsa, by.x="RefSeq", by.y=0, all.x=T)
strain_ann[is.na(strain_ann)] <- 0

## three main cases: (1) empty uniq_strains but not empty blurry_strains; 
## (2) empty blurry but non-empty uniq; (3) non-empty both
## if either of the 4 params are NULL (shared_rrna/prot, unique_rrna/prot) I go to simple collapsing by genus.
## similar logic applies if there are no non-self edges in nodes_and_edges

uniq_strains <- strain_ann[strain_ann$uniq_prot != 0 & strain_ann$prot_genes > 2, ]
blurry_strains <- strain_ann[! strain_ann$RefSeq %in% uniq_strains$RefSeq, ]
strain_ann$score <- strain_ann$uniq_rrna + 10*strain_ann$uniq_prot
merged_clust <- list()

if (is.null(shared_prot) | is.null(shared_rrna) | is.null(uniq_prot) | is.null(uniq_rrna) | nrow(uniq_strains) == 0) { 
  ## case 1: no well defined strains or other gaps, just collapse by genus (each cluster == genus)
  k <- 1
  for (genus in unique(strain_ann$Genus)) { 
    subann <- strain_ann[strain_ann$Genus == genus, ]
    merged_clust[[k]] <- subann$RefSeq
    k <- k + 1
  }
} else { 
  ## case 2: if there is a well-defined network, do the full thing - and it's rather complicated
  edges_all           <- edges_all[,c(1,2,5,6)]
  colnames(edges_all) <- c("node1","node2","overlap","weight")
  edges_all           <- edges_all[edges_all$weight >= 0.5,]
  edges_uniq          <- edges_all[edges_all$node1 %in% uniq_strains$RefSeq & edges_all$node2 %in% uniq_strains$RefSeq, ]
  net_uniq            <- graph_from_data_frame(d=edges_uniq, directed=F)
  net_uniq            <- igraph::simplify(net_uniq, remove.multiple = F) ## this way it does not remove "overlap" from edge attributes!
  lclust              <- cluster_leiden(net_uniq, weights=E(net_uniq)$overlap)

  ### second, we check edge_density of the Leiden clustering and create new (more fine) clusters if there's need
  uniq_clust <- list()
  g2c        <- list()
  
  j <- 1
  for (i in 1:lclust$nb_clusters) { 
    subnet <- induced_subgraph(net_uniq, lclust[[i]])
    dens   <- edge_density(subnet)
    ann    <- strain_ann[strain_ann$RefSeq %in% lclust[[i]], ]
    
    if (dens >= 0.9 | length(lclust[[i]]) == 1) {
      ## if this is a ~clique, we treat it as 1 strain and are happy about it
      uniq_clust[[j]] <- lclust[[i]]
      cat(sprintf("================== cluster %d/%d ::: %d components - %f density===================\n",i,j,nrow(ann),dens))
      print(ann)
      for (genus in unique(ann$Genus)) { 
        g2c[[genus]] <- j
      }
      j <- j + 1
      
    } else {
      ## if this is a crappy cluster, we add top strain per genus
      for (genus in unique(ann$Genus)) {
        subann <- ann[ann$Genus == genus, ]
        uniq_clust[[j]] <- subann$RefSeq
        g2c[[genus]] <- j
        j <- j + 1
        cat(sprintf("================== SUBcluster %d/%d ::: %d components - %f density===================\n",i,j,nrow(subann),dens))
        print(subann)
      }
    }
  }
  
  
  blurry_clust <- list()
  k <- 1
  for (genus in unique(blurry_strains$Genus)) { 
    subann <- blurry_strains[blurry_strains$Genus == genus, ]
    
    ## if we haven's seen this genus yet, make a cluster that's all strains of this genus
    if (is.null(g2c[[genus]])) {
      blurry_clust[[k]] <- subann$RefSeq
      k <- k + 1
      ## if there's a "uniq" cluster with this genus, add the strains to this cluster
    } else { 
      j <- g2c[[genus]]
      uniq_clust[[j]] <- c(uniq_clust[[j]], subann$RefSeq)
    }
  }
  
  ## now the fun part! first let's make a network; edges_all is already filtered to w>0.5
  net_all <- graph_from_data_frame(d=edges_all, directed=F)
  net_all <- igraph::simplify(net_all, remove.multiple = F) ## this way it does not remove "overlap" from edge attributes
  
  ## now let's try to collapse the blurry clusters onto the uniq ones
  ## this happens in 2 cases: 1) if the subgraph made is perfect; 
  ## 2) if the "centroid" (top score) strain matches uniq one very well (>=0.9 overlap from nodes_and_edges.tsv)
  
  b2u <- vector(mode = "list", length(blurry_clust)) ## trick to make a list of NULLs of defined size
  
  for (k in 1:length(blurry_clust)) {
    clique_ovsum <- 0
    clique_clust <- 0 ##finish here
    best_overlap <- 0
    best_weight  <- 0
    best_clust   <- 0
    for (j in 1:length(uniq_clust)) {
      strains <- c(uniq_clust[[j]], blurry_clust[[k]])
      sub1 <- induced_subgraph(net_all, uniq_clust[[j]])
      sub2 <- induced_subgraph(net_all, blurry_clust[[k]]) 
      sub3 <- induced_subgraph(net_all, strains)
      ## count new connections and calculate % of maximum which is n1*n2
      nv1 <- vcount(sub1)
      nv2 <- vcount(sub2)
      ne1 <- ecount(sub1)
      ne2 <- ecount(sub2)
      ne3 <- ecount(sub3)
      ediff <- ne3 - ne1 - ne2
      newconn <- ediff/(nv1*nv2)
      ## calculate sum of all new overlaps
      ovsum <- sum(E(sub3)$overlap) - sum(E(sub2)$overlap) - sum(E(sub1)$overlap)
      
      top1 <- strain_ann[strain_ann$RefSeq %in% blurry_clust[[k]], ]
      top1 <- top1[order(-top1$score), ][1,1]
      top2 <- strain_ann[strain_ann$RefSeq %in% uniq_clust[[j]], ]
      top2 <- top2[order(-top2$score), ][1,1]
      
      edge12 <- try(E(net_all, P = c(top1, top2)), silent=T)
      if (grepl("Error", edge12)) { 
        edge12 <- NULL 
      }
      
      if (newconn >= 0.9 & ovsum > clique_ovsum) {
        cat(sprintf("Clusters: %d blurry, %d uniq, overlap: %d, max: %d\n",k,j,ovsum,clique_ovsum))
        clique_ovsum <- ovsum
        clique_clust <- j
      } else if (! is.null(edge12)) { 
        if (edge12$weight > 0.9 & edge12$overlap > best_overlap) {
          best_overlap <- edge12$overlap
          best_weight <- edge12$weight
          best_clust <- j
        }
      } 
    }
    if (clique_clust != 0) {
      b2u[[k]] <- clique_clust
      cat(sprintf("COLLAPSED as a ~fully connected subnetwork: %d blurry, %d uniq, %d overlap\n",k,clique_clust,clique_ovsum))
    } else if (best_clust != 0) {
      ## second best option
      b2u[[k]] <- best_clust
      cat(sprintf("COLLAPSED by very similar best strains: %d blurry, %d uniq, %f weight, %d overlap\n",k,best_clust,best_weight,best_overlap))
    }
  }
  
  merged_clust <- uniq_clust
  i <- length(uniq_clust) + 1
  for (k in 1:length(blurry_clust)) {
    j <- b2u[[k]]
    if (is.null(j)) { 
      merged_clust[[i]] <- blurry_clust[[k]]
      i <- i + 1
    } else {
      merged_clust[[j]] <- c(merged_clust[[j]], blurry_clust[[k]])
    }
  }
}

r2c <- list()  ## simple refseq ID to final cluster dict
for (i in 1:length(merged_clust)) {
  cat(sprintf("======================= Merged cluster %d ====================\n",i))
  print(strain_ann[strain_ann$RefSeq %in% merged_clust[[i]], ])
  for (refseq in merged_clust[[i]]) {
    r2c[[refseq]] <- i
  }
}

strain_ann$cluster <- as.numeric(lapply(strain_ann$RefSeq, function(x) r2c[[x]]))
strain_ann <- strain_ann[order(strain_ann$cluster, -strain_ann$score, -strain_ann$shared_prot, -strain_ann$shared_rrna), ]
write.table(strain_ann, "filtered.cluster.tsv", sep="\t", quote=F, row.names=F)

