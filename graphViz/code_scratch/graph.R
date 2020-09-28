#figuring out how igraph works 

library(igraph)
library(matie)
library(RANN)
#### SCRATCH ####
admatrix = as.matrix(c(1,2,3))
admatrix = cbind(admatrix,admatrix,admatrix)
admatrix[1,] = c(1,2,3)
admatrix[3,2] = 2

g = graph_from_adjacency_matrix(admatrix,mode = 'directed',weighted = TRUE)
E(g)$width = edge.attributes(g)$weight
plot(g, edge.width = E(g)$width)


plot(layout_with_gem(g))

adj = as.matrix(dist(pcs))
a_sub = matrix(data = 0, nrow = 1368, ncol = 1368)
a_sub[adj>0.45] = adj[adj>0.45]
g1 = graph_from_adjacency_matrix(a_sub,weighted = TRUE)
E(g1)$width = edge.attributes(g1)$weight 
plot(g1,edge.width = E(g)$width)

cs = layout_with_graphopt(g1)
plot(cs)
k = knn(g1)$knn
plot(k)

nn = nn2(t(curr),t(proj),k=1)
nn.idx = nn$nn.idx
nn.dist = nn$nn.dists

cellidx = matrix(data = NA, ncol = 1, nrow = ncells)
gDists = matrix(data = NA, ncol = 1, nrow = ncells)

for (c in seq(1,ncells)){
  print(c)
  query = proj[,c]
  #print(length(query))
  currSpace = curr[,-c]
  #print(dim(currSpace))
  
  currNN = nn2(t(currSpace), query = t(query), k=1)
  #print(currNN)
  #print(currNN$nn.idx)
  #print(currNN$nn.dist)
  
  cellidx[c] = currNN$nn.idx
  gDists[c] = currNN$nn.dist
}




#### SIMULATE ####
obs = t(cbind(c(1,2,3,3,2,1),c(2,1,2,3,4,3)))
exp = t(cbind(c(2,3,3,2,1,1),c(1,2,3,4,3,2)))
par(mfrow = c(1,1))
plot(t(obs),col = "black",main = "Observed",xlim = c(0.5,3.5),ylim = c(0.5,4.5))
text(t(obs)[,1]-0.1,t(obs)[,2],labels = c("A","B","C","D","E","F"), col = "black")
points(t(exp),col = "red",main = "Observed", cex = 2)
text(t(exp)[,1]+0.1,t(exp)[,2],labels = c("A","B","C","D","E","F"), col = "red")
legend(2.5,4.5, legend = c("Observed", "Projected"), pch = c(1,1), col = c("black","red"))
i = sapply(seq(1:ncol(obs)), function(x) nn2(t(obs[,-x]),t(exp[,x]),k=1)$nn.idx)
d = sapply(seq(1:ncol(obs)), function(x) nn2(t(obs[,-x]),t(exp[,x]),k=1)$nn.dist)

w = 1/(1+d)

for (c in seq(1,length(i))){
  if (i[c]>=c){
    i[c] = i[c] + 1
  }
}


el = cbind(seq(1,ncol(obs)),i)
gsim = graph_from_edgelist(el,directed =TRUE)
gsimFD = layout_with_fr(gsim)
par(mfrow = c(1,2))
plot(gsim, main = "Default directed graph")
plot(gsimFD, main = "FDG")


#### U2O5 DATA ####
rvel.cd = readRDS("rvelcd.rds")

curr = rvel.cd$current
proj = rvel.cd$projected


#find nearest neighbor to each projected cell in current cell 
#build a graph where edges are pointing from the cell to the nearest neighbor of its projected cell 
#need to exclude current cell 

ncells = ncol(curr)

cellidx = sapply(seq(1:ncells), function(x) nn2(t(curr[,-x]),t(proj[,x]),k=1)$nn.idx)
celldist = sapply(seq(1:ncells), function(x) nn2(t(curr[,-x]),t(proj[,x]),k=1)$nn.dist)


for (c in seq(1,length(cellidx))){
  if (cellidx[c]>=c){
    cellidx[c] = cellidx[c] + 1
  }
}

weights = 1/(1+gDists)

par(mfrow=c(1,1))
edgeList = cbind(seq(1,ncells),cellidx)
g = graph_from_edgelist(edgeList,directed =TRUE)
plot(g)
edge.attributes(g)$weight = weights
E(g)$width = edge.attributes(g)$weight
plot(g, edge.width = E(g)$width)


gFD = layout_with_fr(g)
plot(gFD)

colnames(gFD) = c("C1","C2")
rownames(gFD) = colnames(curr)


show.velocity.on.embedding.cor(scale(gFD), rvel.cd, n=100, scale='sqrt', cell.colors=cell.color,
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2)


#starting coordinates = PCA  -- doesn't work 
gFD.pca = layout_with_fr(g,emb.test)
plot(gFD.pca)

colnames(gFD.pca) = c("C1","C2")
rownames(gFD.pca) = colnames(curr)

show.velocity.on.embedding.cor(scale(gFD.pca), rvel.cd, n=100, scale='sqrt', cell.colors=cell.color,
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2)
