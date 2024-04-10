##library(ape)
##library(Matrix)
library(phybase)
sptree <- "((((((((((M3:1.0,M4:1.0):1.0,(M1:1.0,M2:1.0):1.0):1.0,M5:3.0):1.0,M6:4.0):1.0,(((M8:1.0,M9:1.0):1.0,M7:2.0):1.0,M10:3.0):2.0):1.0,M11:6.0):1.0,M12:7.0):1.0,(M13:1.0,M14:1.0):7.0):1.0,(((M17:1.0,M16:1.0):1.0,M15:2.0):1.0,M18:3.0):6.0):1.0,OUTGROUP:10.0);"  #不需要支持度,且需等枝长
spname <- species.name(sptree)
nodematrix <- read.tree.nodes(sptree,spname)$nodes
nodematrix[,5] <- 20  #需测试出合适的值，一般以最小枝长的4倍及以上为宜
#node.height(36,nodematrix,19) #36为物种数+节点数-1，19为物种数
rootnode <- rootoftree(nodematrix)
nspecies <- length(spname)
seq <- rep(1,nspecies)
gentrees <- rep("",10000)
for(i in 1:10000){
  gentrees[i] <- sim.coaltree.sp(rootnode,nodematrix,nspecies,seq,name=spname)$gt
  }
  write(gentrees,"simulated.tre")
