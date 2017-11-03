library(ggtree)
nwk<- "/Users/ryannguyen/Desktop/Math_127/Phylogenetics/tree.nwk"
tree<- read.tree(nwk)
ggtree(tree, color="firebrick", size=1) + geom_tiplab() +geom_tippoint()

