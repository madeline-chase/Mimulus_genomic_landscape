### This script is used to calculate the tree concordance scores between the genome-wide tree and the genomic windows trees, by calculating the correlation coefficient between the distance matrices.



library(ape)
library(MASS)
## Load genome-wide tree and convert to distance matrix
species_tree <- read.tree('trees/RAxML_bestTree.WG_tree_G1_9_taxa')
species_dist_mat <- cophenetic.phylo(species_tree)
species_dist_ordered <- species_dist_mat[order(rownames(species_dist_mat)),]
species_dist_ordered <- species_dist_ordered[,order(colnames(species_dist_ordered))]

## Load window trees
win_trees <- read.tree('trees/window_trees/win_tree_500kb_tops_only_ends_rmvd.txt')
dist_cors_list <- c()


## Get correlation between window tree and genome-wide tree for each window
for(i in 1:length(win_trees)){
  dist_mat_win <- cophenetic.phylo(win_trees[[i]])
  dist_mat_ordered <- dist_mat_win[order(rownames(dist_mat_win)),]
  dist_mat_ordered <- dist_mat_ordered[,order(colnames(dist_mat_ordered))]
  dist_mat_species_cor <- cor(c(dist_mat_ordered), c(species_dist_ordered))
  dist_cors_list[[i]] <- dist_mat_species_cor
}

## Write concordance scores to file
dist_cors_mat <- as.matrix(dist_cors_list)
write.matrix(dist_cors_mat, file='window2species_corrs_500kb.txt', sep = '\t')



