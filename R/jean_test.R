## setup
library(Rcpp)
library(reticulate)
library(igraph)
library(matie)
library(RANN)
library(velocyto.R)
library(MUDAN)
library(RSpectra)
library(umap)
source("../graphViz/projectedNeighbors.R")
sourceCpp("../graphViz/graphVizC.cpp")

## data
use_condaenv("velo",required = TRUE)
scv = import("scvelo")
adata = scv$datasets$dentategyrus()
#adata = scv$datasets$pancreas()

## old embedding
emb_umap = adata$obsm['X_umap'] #extract umap embedding
clusters = adata$obs$clusters #extract clusters
rownames(emb_umap) = names(clusters) = adata$obs_names$values
col = rainbow(length(levels(clusters)),s = 0.8, v = 0.8)
cell.cols = col[clusters] #color according to cluster
names(cell.cols) = names(clusters)
plotEmbedding(emb_umap, groups=clusters, xlab = "UMAP X", ylab = "UMAP Y", mark.clusters=TRUE)

## get counts
spliced = as.matrix(t(adata$layers['spliced']))
unspliced = as.matrix(t(adata$layers['unspliced']))
cells = adata$obs_names$values
genes = adata$var_names$values
colnames(spliced) = colnames(unspliced) = cells
rownames(spliced) = rownames(unspliced) = genes


##### even out clusters 
table(clusters)
big.clusters = c("Astrocytes","Granule immature","Granule mature", "Neuroblast") #>100 cells 
clusters.even = clusters[!(clusters %in% big.clusters)]
for (b in big.clusters){
  curr.big.cluster = clusters[clusters==b]
  currN = length(curr.big.cluster)
  sub  = curr.big.cluster[seq(1,currN,currN/50)]
  clusters.even = c(clusters.even,sub)
}

spliced = spliced[,names(clusters.even)]
unspliced = unspliced[,names(clusters.even)]
cells = colnames(spliced)
cell.cols = cell.cols[cells]

#####



## filter
gexpS = log10(rowSums(spliced)+1)
gexpU = log10(rowSums(unspliced)+1)

par(mfrow = c(1,2))
hist(gexpS, breaks = 500, main = "log spliced counts")
hist(gexpU, breaks = 500, main = "log unspliced counts")

## less stringent?
goodGenes = genes[gexpS > 0 & gexpU > 0]
spliced = spliced[goodGenes,]
unspliced = unspliced[goodGenes,]

dim(spliced)

## normalize
allExp = spliced + unspliced #use combined spliced and unspliced counts
all.cpm = MUDAN::normalizeCounts(allExp) #cpm normalize
all.norm = MUDAN::normalizeVariance(all.cpm,details = TRUE, plot = FALSE) #variance stabilize
odsGenes = goodGenes[all.norm$ods] #overdispersed genes
all.logODS = log10(as.matrix(t(all.norm$mat[all.norm$ods,]))+1) #keep overdispersed genes and log normalize

## PCA
pca = svds(A = all.logODS, k=50, opts = list(center = TRUE, scale = FALSE, maxitr = 2000, tol = 1e-10))
var = pca$d
plot(var)

## plot
emb.pca = pca$u[,1:2]
row.names(emb.pca) = row.names(all.logODS)
plotEmbedding(emb.pca, groups=clusters, xlab = "PC1", ylab = "PC2", main = "PCA on subsampled data", mark.clusters = TRUE)

#calculate cell-cell distance in pc space for velocity calculation
cell.dist = as.dist(1-cor(t(pca$u))) #cell distance in PC space

## umap
#pcsToKeep = pca$u[,1:15]
pcsToKeep = pca$u[,1:30]
set.seed(1)
library(uwot)  
emb.umap = umap(pcsToKeep, min=1) ## try to make it look closer to example
row.names(emb.umap) = row.names(all.logODS)
par(mfrow=c(1,1))
plotEmbedding(scale(emb.umap), groups=clusters, xlab = "UMAP X", ylab = "UMAP Y",main = "UMAP", mark.clusters = TRUE)
plotEmbedding(scale(emb_umap), groups=clusters, xlab = "UMAP X", ylab = "UMAP Y",main = "UMAP scvelo", mark.clusters = TRUE)

## tSNE
set.seed(1)
library(Rtsne)
emb.tsne = Rtsne(pcsToKeep, k=10)$Y
row.names(emb.tsne) = row.names(all.logODS)
plotEmbedding(scale(emb.tsne), groups=clusters, xlab = "tSNE X", ylab = "tSNE Y",main = "tSNE", mark.clusters = TRUE)

## diffusion maps
library(destiny)
emb.dm <- DiffusionMap(pcsToKeep)
emb.dm <- eigenvectors(emb.dm)[,1:2]
row.names(emb.dm) = row.names(all.logODS)
plotEmbedding(scale(emb.dm), groups=clusters, xlab = "DC X", ylab = "DC Y",main = "Diffusion Map", mark.clusters = TRUE)





## load model
vel = readRDS("../graphViz/saved_objects/all_cells_models/neuro_vel_k30.rds")
#vel = readRDS("../graphViz/saved_objects/all_cells_models/panc_vel_k30.rds")
vel$kCells
curr = vel$current
proj = vel$projected

##### even out clusters 

curr = curr[,names(clusters.even)]
proj = proj[,names(clusters.even)]



##### 



## show velocity
# show.velocity.on.embedding.cor(scale(emb.pca), vel, n=100,
#                                scale='sqrt', cell.colors=cell.cols,
#                                cex=1, arrow.scale=1, show.grid.flow=TRUE,
#                                min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2,
#                                main = "Velocities on PCA embedding")
# 
# show.velocity.on.embedding.cor(scale(emb.umap), vel, n=100,
#                                scale='sqrt', cell.colors=cell.cols,
#                                cex=1, arrow.scale=2, show.grid.flow=TRUE,
#                                min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2,
#                                main = "Velocities on UMAP embedding")
# 
# show.velocity.on.embedding.cor(scale(emb_umap), vel, n=100,
#                                scale='sqrt', cell.colors=cell.cols,
#                                cex=1, arrow.scale=2, show.grid.flow=TRUE,
#                                min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2,
#                                main = "Velocities on original UMAP embedding")
# 
# show.velocity.on.embedding.cor(scale(emb.tsne), vel, n=100,
#                                scale='sqrt', cell.colors=cell.cols,
#                                cex=1, arrow.scale=2, show.grid.flow=TRUE,
#                                min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2,
#                                main = "Velocities on tsne embedding")


### try graph from just pcs?
ks = c(5,10,20,50)
pcs <- pcsToKeep
par(mfrow = c(2,2))
for (k in ks){
  nn = RANN::nn2(pcs, k = k) ## KNN
  names(nn) <- c('idx', 'dists')
  weight <- 1/(1+ as.vector(nn$dists))
  nn.df = data.frame(from = rep(1:nrow(nn$idx), k),
                     to = as.vector(nn$idx),
                     weight = weight
  )
  g <- igraph::graph_from_data_frame(nn.df, directed = FALSE)
  g <- igraph::simplify(g)
  fdg = layout_with_fr(g,dim=2)
  colnames(fdg) = c("C1","C2")
  rownames(fdg) = row.names(all.logODS)
  plotEmbedding(scale(fdg), groups=clusters, xlab = "fdg X", ylab = "fdg Y",main = paste("fdg k:",k), mark.clusters = TRUE)
  
}


# show.velocity.on.embedding.cor(scale(fdg), vel, n=100,
#                                scale='sqrt', cell.colors=cell.cols,
#                                cex=1, arrow.scale=2, show.grid.flow=TRUE,
#                                min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2,
#                                main = "Velocities on fdg embedding", n.cores=10)

##### velocity-informed fdg
u = pca$u #scores
v = pca$v #loads
row.names(v) = odsGenes

common.genes = sapply(row.names(curr), function(x) x %in% odsGenes)
common.genes.names = rownames(curr)[common.genes]
curr.pca = curr[common.genes,]
v.sub = v[row.names(curr.pca),]

par(mfrow = c(1,2))

curr.scores = t(curr.pca) %*% v.sub
plot(curr.scores[,1:2], col = cell.cols, pch = 16, main = 'Observed')
legend(x=2, y=1.2, legend = levels(clusters), col = col, pch =16, cex = 0.8)

proj.pca = proj[common.genes,]
proj.scores = t(proj.pca) %*% v.sub
plot(proj.scores[,1:2], col = cell.cols, pch = 16, main = 'Projected') #,xlim = c(-0.5,4),ylim = c(-3,1.25)
legend(x=2, y=1.2, legend = levels(clusters), col = col, pch =16, cex = 0.8)

set.seed(1)
cell.names = colnames(curr)
#cell.names.sub = cell.names[sample(c(1:length(cell.names)),400)] #downsample cells for speed
cell.names.sub = cell.names #no downsampling
ncells.sub = length(cell.names.sub)
cell.cols.grph = cell.cols[cell.names.sub]

curr.scores.cellsub = t(curr.scores[cell.names.sub,]) #### change here to change number of PCs included
proj.scores.cellsub = t(proj.scores[cell.names.sub,]) ####
#plot(t(curr.scores.cellsub[1:2,]), col = cell.cols.grph, pch = 16)

par(mfrow = c(2,2))
ks = c(5,10,20,50)
for (k in ks){
  t = 0.4
  set.seed(1)
  velograph = graphViz(curr.scores.cellsub,
                       proj.scores.cellsub,
                       k,"L2","cosine",t,
                       weighted = TRUE,
                       cell.colors = cell.cols.grph,
                       title = paste("K =",k,"weighted, thresh =",t),
                       plot=TRUE, return_graph=TRUE)
  
  #plotEmbedding(velograph$fdg_coords, groups=clusters, xlab = "fdg X", ylab = "fdg Y",main = paste("fdg k:",k), alpha = 0.8, mark.clusters = TRUE, mark.cluster.cex = 1)
}



#compare pca graph to velocity graph
par(mfrow = c(1,2))
plotEmbedding(scale(fdg), groups=clusters, xlab = "fdg X", ylab = "fdg Y",main = "fdg", mark.clusters = TRUE)
plotEmbedding(velograph$fdg_coords, groups=clusters, xlab = "fdg X", ylab = "fdg Y", mark.clusters=TRUE)


## compare simple weights versus velocity informed weights
par(mfrow = c(1,1))
plot(weight, E(velograph$graph)$weight, pch = 4, xlab = "simple weights", ylab = "velocity weights")

## even out clusters 









##### Scratch
par(mfrow = c(2,1))
hist(velograph$projected_neighbors$edge_weights, breaks = 500)
hist(E(velograph$graph)$weight, breaks = 500)
par(mfrow = c(1,1))
plot(velograph$projected_neighbors$edge_weights,E(velograph$graph)$weight)
hist(velograph$projected_neighbors$dist_comp$negVectSim,breaks = 500)
hist(log(velograph$projected_neighbors$dist_comp$invDist) ,breaks = 500)
hist(velograph$projected_neighbors$all_dists, breaks = 500, ylim = c(0,1e01))
hist(velograph$projected_neighbors$edge_weights, breaks = 500)
a = velograph$projected_neighbors$edge_weights
hist(a, breaks = 500)
b = max(a) - a
hist(b, breaks = 500)

show.velocity.on.embedding.cor(scale(velograph$fdg_coords), vel, n=100,
                               scale='sqrt', cell.colors=cell.cols,
                               cex=1, arrow.scale=2, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2,
                               main = "Velocities on velograph embedding", n.cores=10)









################## Building on Lyla's simple simulation
set.seed(1)
obs = jitter(t(cbind(c(1,2,3,3,2,1, 2, 2, 2),c(2,1,2,3,4,3,0, -1, -2))))
exp = jitter(t(cbind(c(2,3,3,2,1,1, 2, 2, 2),c(1,2,3,4,3,2,1,0,-1))))
labels = c("A","B","C","D","E","F","X","Y","Z")
par(mfrow = c(1,1))
plot(t(obs),col = "black",main = "Sim", pch = 16)
text(t(obs)[,1]-0.1,t(obs)[,2],labels = labels, col = "black")
points(t(exp),col = "red", cex = 2)
text(t(exp)[,1]+0.1,t(exp)[,2],labels = labels, col = "red")
legend(2.5,4.5, legend = c("Observed", "Projected"), pch = c(1,1), col = c("black","red"))

i = sapply(seq(1:ncol(obs)), function(x) nn2(t(obs[,-x]),t(exp[,x]),k=1)$nn.idx)
d = sapply(seq(1:ncol(obs)), function(x) nn2(t(obs[,-x]),t(exp[,x]),k=1)$nn.dist)
w = 1/(1+d)

#since current observed cell is being excluded when searching for nearest neighbor
#any index >= i will be off by 1.
for (c in seq(1,length(i))){
  #correcting indices
  if (i[c]>=c){
    i[c] = i[c] + 1
  }
}

el = cbind(seq(1,ncol(obs)),i) #edge list
rownames(el) <- labels
gsim = graph_from_edgelist(el,directed =TRUE)
V(gsim)$label <- labels
V(gsim)$size = 10

set.seed(1)
gsimFD = layout_with_fr(gsim)
par(mfrow = c(1,2))
plot(gsim, main = "Default directed graph")
plot(gsimFD, main = "FDG: vertex coordinates")
text(gsimFD+0.1, labels = labels)

library(umap)
test <- umap(el, n_neighbors = 2)
plot(test, main = "umap")
text(test+0.1, labels=labels)

test <- prcomp(el, scale=TRUE, center=TRUE)
test <- test$x
plot(test, main = "pca")
text(test+0.1, labels=labels)

test <- Rtsne(el, perplexity=2)$Y
plot(test, main = "tsne")
text(test+0.1, labels=labels)

nn = RANN::nn2(el, k = 2) ## KNN
names(nn) <- c('idx', 'dists')
weight <- 1/(1+ as.vector(nn$dists))
nn.df = data.frame(from = rep(1:nrow(nn$idx), k),
                   to = as.vector(nn$idx),
                   weight = weight
)
g <- igraph::graph_from_data_frame(nn.df, directed = FALSE)
g <- igraph::simplify(g)
fdg = layout_with_fr(g,dim=2)
colnames(fdg) = c("C1","C2")
rownames(fdg) <- labels
plot(fdg, main = "simple fdg")
text(fdg+0.1, labels = labels)
