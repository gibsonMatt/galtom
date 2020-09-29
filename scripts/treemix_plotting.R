#setwd("~/Box\ Sync/Projects/Galapagos/Sequencing/pipeline2/dataset1/scripts/treemix")

###################Imports###########################
#####################################################

library(ape)
library(phangorn)
library(phytools)

###################Data##############################
#####################################################

consensus <- read.tree(text="(MG115:0.0001,((MG120_G:0.00859751,MG120_C:0.00626863):0.0122869,(((((((MG114:0.00100054,MG105:0.00480819):0.00133275,MG117:0.0001):0.00098446,MG113:0.00204414):0.000718844,MG111:0.0122304):0.0047877,MG107:0.00501234):0.0303889,ecu:0.0001):0.00516695,peru:0.0537063):0.0488039):0.0001);")
phy <- read.tree(file = "bootstrap_trees1.txt") #bootstrapped trees. 500 SNPs sampled with replacement. 1,000 replicates

#####################################################
#####################################################

#smooth the tree with chronopl
#pdf("smoothed_consensus.pdf", 5,6)
plot(chronopl(consensus, lambda = 1))
#dev.off()

smoothed <- list(consensus)
class(smoothed) <- "multiPhylo"

#Smooth all the bootstrapped trees
for (y in phy){
  #x <- reroot(y, 15, 0.01)
  x <- chronopl(y, lambda = 1)
  smoothed <- c(smoothed, x)
}


#####################################################
#####################################################


#DensiTree plot
par(mfrow=c(1,1))
#pdf("smoothed_reps.pdf", 5,6)
densiTree(smoothed, alpha=0.05, col = "black", consensus = consensus, scaleX = F, add=TRUE)
#dev.off()

#####################################################
#####################################################


is.monophyletic(consensus, c('MG115', 'MG120_G', 'MG120_C', "ecu"), plot=T)

#x <- is.monophyletic(p, c('MG107', 'MG111', 'MG113', 'MG117', 'MG105', 'MG114'), plot=F)

count <- 0
tote <- 0
for (p in smoothed){
  tote <- tote+1
  x <- is.monophyletic(p, c('ecu', 'peru'), plot=F)
  if (x == TRUE){
    count <- count + 1
  }
}
count/tote

count <- 0

for (p in smoothed){
  x <- is.monophyletic(p, c('MG116', 'MG117', 'MG114', 'MG113'), plot=F)
  if (x == TRUE){
    count <- count + 1
  }
}
count/1000

un <- unique.multiPhylo(smoothed)

#write.tree(consensus, file = "../ml_tree/consensus_tree.tree")
#write.tree(smoothed, file = "../ml_tree/gene_trees.tree")
