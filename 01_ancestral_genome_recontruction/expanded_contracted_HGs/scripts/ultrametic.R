# Load necessary libraries
install.packages("ape")
install.packages("TreeTools")
install.packages("phytools")

library(ape)        
library(TreeTools) 
library(phytools)  


mytree <- read.tree("Terrestrialization/phylogeny_tree/species_tree_build/iqtree_results/iqtree_C60_results_ConstrainTree/C60_CT_metazoa.treefile")

# Root the tree
rooted <- RootTree(mytree, c('Aque', 'Emue'))

# Convert the tree to binary
binTree <- MakeTreeBinary(rooted)

# Convert the binary tree to ultrametric
ultrametric <- force.ultrametric(binTree, method = "extend")

# Save the ultrametric tree in Newick format
write.tree(ultrametric, "Terrestrialization/expansion/cafe/ultrametric_metazoa_tree.newick")

# Optional: Plot the ultrametric tree
plot(ultrametric)
