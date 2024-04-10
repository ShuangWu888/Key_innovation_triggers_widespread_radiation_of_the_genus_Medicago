~/software/julia-1.6.2/bin/julia -p 10
using Pkg
using PhyloNetworks
using PhyloPlots
using RCall

iqtreeCF=readTrees2CF("2.arstral.CDS.simple.reroot.order.rename.tree")
tre=readTopology("4.arstral.CDS.simple.tree.BS10.individual.reroot.order.x.order.newick")

net0 = snaq!(tre,  iqtreeCF, hmax=0, filename="snaq/net0_iqtreeCF", seed=123, runs=100)
net1 = snaq!(net0, iqtreeCF, hmax=1, filename="snaq/net1_iqtreeCF", seed=123, runs=100)
net2 = snaq!(net1, iqtreeCF, hmax=2, filename="snaq/net2_iqtreeCF", seed=123, runs=100)
net3 = snaq!(net2, iqtreeCF, hmax=3, filename="snaq/net3_iqtreeCF", seed=123, runs=100)
