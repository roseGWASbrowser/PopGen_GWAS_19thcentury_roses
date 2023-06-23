# TL - 010222
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/first_round_calling/kinship")
IDs=read.table("IDs4kinshipnetwork",h=T)
relationships=read.table("relationships4network.txt",h=T)

####
require(igraph)

par(mar=c(0,0,0,0))
g=graph.data.frame(relationships, directed=FALSE)
graph_layout=layout.fruchterman.reingold(g)
plot(g,layout=graph_layout)

plot(g,layout=graph_layout,edge.width=4-relationships$Weight)  # ,vertex.color=IDs$Type


net = graph_from_data_frame(relationships,vertices = IDs,directed = F)
colrs.v = c(AncientAsia = "red", AncientEU = "gold", Botanical = "forestgreen", FirstHyb = "lightblue", HybTea = "pink",Ambiguous="white") #node colours
V(net)$color = colrs.v[V(net)$Type]

colrs.e = c(level1 = "firebrick1", level2 = "darkorange",level3="gold1") #edge colours
E(net)$color = colrs.e[E(net)$Weight] 

plot(net, edge.curved=0.2,vertex.size=10,edge.width=3)


