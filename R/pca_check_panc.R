########### Get Pancreas data
library(reticulate)
conda_list()
use_condaenv("r-velocity", required = TRUE)
scv = import("scvelo")

adata = scv$datasets$pancreas()

emb_umap = adata$obsm['X_umap'] #extract umap embedding
clusters = adata$obs$clusters #extract clusters
rownames(emb_umap) = names(clusters) = adata$obs_names$values

col = rainbow(length(levels(clusters)),s = 0.8, v = 0.8)
cell.cols = col[clusters] #color according to cluster
names(cell.cols) = names(clusters)

plot(emb_umap, col = cell.cols, pch=16, xlab = "UMAP X", ylab = "UMAP Y", xlim = c(-15,20))
legend(x=10, y=6, legend = levels(clusters), col = col, pch = 16)

spliced = as.matrix(t(adata$layers['spliced']))
unspliced = as.matrix(t(adata$layers['unspliced']))
cells = adata$obs_names$values
genes = adata$var_names$values
colnames(spliced) = colnames(unspliced) = cells
rownames(spliced) = rownames(unspliced) = genes

#### velocity
#vel <- gene.relative.velocity.estimates(spliced, unspliced,
#                                        deltaT=1, kCells=30,
#                                        cell.dist=cell.dist)
vel = readRDS('../panc_vel_k30.rds')
vel$kCells
curr = vel$current
proj = vel$projected

## normalize variance on current
## use same normalizing factor for projected
foo <- normalizeVariance(curr, details=TRUE, plot=TRUE, alpha = 0.2) ## more relaxed threshold
names(foo)
head(foo$df)
## double check
ods.genes <- rownames(foo$mat)[foo$ods]

par(mfrow=c(2,2))
plot(foo$df$m, foo$df$v)
points(foo$df[ods.genes,]$m, foo$df[ods.genes,]$v, col='red')
plot(foo$df$m, foo$df$res)
points(foo$df[ods.genes,]$m, foo$df[ods.genes,]$res, col='red')

## variance normalize both
## the current and projected transcriptional states
## by same scaling factor
curr.sub <- curr[ods.genes,] * foo$df[ods.genes,]$gsf
table(curr.sub[,1] == foo$mat[ods.genes,1]) ## double check
proj.sub <- proj[ods.genes,] * foo$df[ods.genes,]$gsf

plot(apply(curr.sub, 1, mean), apply(curr.sub, 1, var))
plot(apply(proj.sub, 1, mean), apply(proj.sub, 1, var))

curr.sub <- log10(curr.sub+1)
proj.sub <- log10(proj.sub+1)

library(RSpectra)
pca = svds(A = t(curr.sub), k=50, opts = list(center = TRUE, scale = FALSE, maxitr = 2000, tol = 1e-10))
range(pca$u)
var = pca$d
plot(var)

nPCs <- 10

m <- t(curr.sub)
m <- m - apply(m, 1, mean)
range(apply(m, 1, mean))
pca.curr <- m %*% pca$v[,1:nPCs]
rownames(pca.curr) <- colnames(curr.sub)
head(pca.curr)

m <- t(proj.sub)
m <- m - apply(m, 1, mean)
range(apply(m, 1, mean))
pca.proj <- m %*% pca$v[,1:nPCs]
rownames(pca.proj) <- colnames(proj.sub)
head(pca.proj)

colnames(pca.curr) <- colnames(pca.proj) <- paste0('PC', 1:ncol(pca.curr))

## double check zero centered
par(mfrow=c(2,2))
hist(pca.curr)
hist(pca.proj)

################ our visualization
library(igraph)
library(matie)
library(RANN)
library(Rcpp)
getwd()
setwd('graphViz')
#source('projectedNeighbors.R')
source("projectedNeighbors_weightedCD.R")

gsim.gg = graphViz(t(pca.curr), t(pca.proj), 30,
                   cell.colors=cell.cols,
                   similarity_threshold=0.1,
                   distance_weight = 2,
                   weighted=FALSE,
                   plot = FALSE,
                   return_graph = TRUE)
par(mfrow=c(1,1))
plot(gsim.gg$fdg_coords, main = "FDG: vertex coordinates", col=cell.cols, pch=16)

par(mfrow=c(1,1))
g <- gsim.gg$graph
V(g)$label = NA
V(g)$size = 2
E(g)$arrow.size = 1
E(g)$color = rgb(0,0,0,0.01)
plot(g)

deltaExp = vel$deltaE
scores = consistency(gsim.gg$fdg_coords, deltaExp, 30, FALSE)
mean(scores)

############### Explore parameter
results <- do.call(rbind, lapply(c(1,2,3,5,10,15,30,50), function(k) {
  print(k)
  sapply(seq(0,2,0.25), function(w) {
    print(w)
    gsim.gg = graphViz(t(pca.curr), t(pca.proj), k,
                       cell.colors=cell.cols,
                       similarity_threshold=0.1,
                       distance_weight = w,
                       weighted=FALSE,
                       plot = FALSE,
                       return_graph = TRUE)
    deltaExp = vel$deltaE
    scores = consistency(gsim.gg$fdg_coords, deltaExp, 30, FALSE)
    mean(scores)
  })
}))

range(results)
results[1:5,1:5]
rownames(results) <- paste0('K:', c(1,2,3,5,10,15,30,50))
colnames(results) <- paste0('w:', seq(0,2,0.25))

heatmap(results, Rowv=NA, Colv=NA, scale='none', col=colorRampPalette(c('white', 'black'))(100))

heatmap(results[4:8,2:9], Rowv=NA, Colv=NA, col=colorRampPalette(c('white', 'black'))(100), scale='none')
heatmap(results[4:8,2:9], Rowv=NA, Colv=NA, col=colorRampPalette(c('white', 'black'))(100), scale='row')
heatmap(results[4:8,2:9], Rowv=NA, Colv=NA, col=colorRampPalette(c('white', 'black'))(100), scale='col')
range(results[4:8,2:9])

## which is best
which(results == max(results), arr.ind=TRUE)
rownames(results)[5]
colnames(results)[6]

## are we that much worse off if we don't consider weight
which(results == max(results[, 'w:1']), arr.ind=TRUE)
rownames(results)[8]
colnames(results)[5]

gsim.gg = graphViz(t(pca.curr), t(pca.proj), 50,
                   cell.colors=cell.cols,
                   similarity_threshold=0.1,
                   distance_weight = 1,
                   weighted=FALSE,
                   plot = FALSE,
                   return_graph = TRUE)
par(mfrow=c(1,1))
plot(gsim.gg$fdg_coords, main = "FDG: vertex coordinates", col=cell.cols, pch=16)

par(mfrow=c(1,1))
g <- gsim.gg$graph
V(g)$label = NA
V(g)$size = 2
E(g)$arrow.size = 1
E(g)$color = rgb(0,0,0,0.01)
plot(g)

deltaExp = vel$deltaE
scores = consistency(gsim.gg$fdg_coords, deltaExp, 30, FALSE)
## global mean
mean(scores)

## mean per cluster
foo <- model.matrix(~ 0 + clusters)
bar <- t(foo) %*% scores
scores.cluster <- bar/colSums(foo)
scores.cluster

###### Assess stability across clusters as function of K
results2 <- do.call(rbind, lapply(c(1,2,3,5,10,15,30,50), function(k) {
  print(k)
    gsim.gg = graphViz(t(pca.curr), t(pca.proj), k,
                       cell.colors=cell.cols,
                       similarity_threshold=0.1,
                       distance_weight = 1,
                       weighted=FALSE,
                       plot = FALSE,
                       return_graph = TRUE)
    deltaExp = vel$deltaE
    scores = consistency(gsim.gg$fdg_coords, deltaExp, 30, FALSE)
    foo <- model.matrix(~ 0 + clusters)
    bar <- t(foo) %*% scores
    scores.cluster <- bar/colSums(foo)
    t(scores.cluster)
}))

heatmap(results2, Rowv=NA, Colv=NA, scale='none', col=colorRampPalette(c('white', 'black'))(100))





