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
net4 = snaq!(net3, iqtreeCF, hmax=4, filename="snaq/net4_iqtreeCF", seed=123, runs=100)
net5 = snaq!(net4, iqtreeCF, hmax=5, filename="snaq/net5_iqtreeCF", seed=123, runs=100)
#net6 = snaq!(net5, iqtreeCF, hmax=6, filename="snaq/net6_iqtreeCF", seed=123, runs=100)
#net7 = snaq!(net6, iqtreeCF, hmax=7, filename="snaq/net7_iqtreeCF", seed=123, runs=100)
#net8 = snaq!(net7, iqtreeCF, hmax=8, filename="snaq/net8_iqtreeCF", seed=123, runs=100)
#net9 = snaq!(net8, iqtreeCF, hmax=9, filename="snaq/net9_iqtreeCF", seed=123, runs=100)
#net10 = snaq!(net9, iqtreeCF, hmax=10, filename="snaq/net10_iqtreeCF", seed=123, runs=100)

R"pdf"("plot-all-net.pdf")  #, width=3, height=3);
#plot(net0, :R);
R"layout(matrix(c(1,2,3,4,5,6,7,8,0),3,3,byrow = T))"
R"par"(mar=[0.1,0.1,0.1,0.5])
#rootatnode!(net0,10)
plot(net0, :R,showGamma=true)
plot(net1, :R,showGamma=true)
plot(net2, :R,showGamma=true)
plot(net3, :R,showGamma=true)
plot(net4, :R,showGamma=true)
plot(net5, :R,showGamma=true)
#plot(net6, :R,showGamma=true)
#plot(net7, :R,showGamma=true)
#plot(net8, :R,showGamma=true)
#plot(net9, :R,showGamma=true)
#plot(net10, :R,showGamma=true)

R"dev.off()";
