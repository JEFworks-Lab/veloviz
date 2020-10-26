##### Overdispersed genes from full dataset, not just velocity-informative genes

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
names(vel)

vel$kCells
curr = vel$current
proj = vel$projected

mat <- vel$conv.emat.norm + vel$conv.nmat.norm
range(curr)
range(proj)
range(mat)

## normalize variance on current
## use same normalizing factor for projected
foo <- normalizeVariance(mat, details=TRUE, plot=TRUE, alpha = 0.05)
names(foo)
head(foo$df)
## double check
ods.genes <- rownames(foo$mat)[foo$ods]
ods.genes <- ods.genes2 <- intersect(ods.genes, rownames(curr)) ## limit to velocity informative genes
length(ods.genes)

## regular analysis
test <- log10(foo$mat[ods.genes, ]+1)
pca = svds(A = t(test), k=50, opts = list(center = TRUE, scale = TRUE, maxitr = 2000, tol = 1e-10))
#pca = svds(A = t(test), k=50, opts = list(center = TRUE, scale = FALSE, maxitr = 2000, tol = 1e-10))

m <- t(test)
m <- m - apply(m, 1, mean)
range(apply(m, 1, mean))
m <- m / apply(m, 1, sd)
range(apply(m, 1, var))
nPCs <- 10
pca.test <- m %*% pca$v[,1:nPCs]
emb_test <- umap(pca.test, min_dist = 0.5)
plot(emb_test, main = "Umap test", col=cell.cols, pch=16)




# ## alternative variance normalize on curr
# foo <- normalizeVariance(curr, details=TRUE, plot=TRUE, alpha = 0.05)
# names(foo)
# head(foo$df)
# ## double check
# ods.genes <- rownames(foo$mat)[foo$ods]
# length(ods.genes)
# ods.genes <- ods.genes2 <- intersect(ods.genes, rownames(curr)) ## limit to velocity informative genes
# length(ods.genes)
#
# ## establish PCs on curr too
# test <- log10(foo$mat[ods.genes, ]+1)
# #pca = svds(A = t(test), k=50, opts = list(center = TRUE, scale = TRUE, maxitr = 2000, tol = 1e-10))
# pca = svds(A = t(test), k=50, opts = list(center = TRUE, scale = FALSE, maxitr = 2000, tol = 1e-10))
#
# nPCs <- 10
# m <- t(test)
# m <- m - apply(m, 1, mean)
# range(apply(m, 1, mean))
# #m <- m / apply(m, 1, sd)
# #range(apply(m, 1, var))
# pca.test <- m %*% pca$v[,1:nPCs]
# emb_test <- umap(pca.test, min_dist = 0.5)
# plot(emb_test, main = "Umap test", col=cell.cols, pch=16)


## variance normalize both
## the current and projected transcriptional states
## by same scaling factor
#ods.genes2 <- intersect(ods.genes, rownames(curr)) ## limit to velocity informative genes
curr.sub <- curr[ods.genes2,] * foo$df[ods.genes2,]$gsf
proj.sub <- proj[ods.genes2,] * foo$df[ods.genes2,]$gsf

par(mfrow=c(2,2))
plot(apply(curr.sub, 1, mean), apply(curr.sub, 1, var))
plot(apply(proj.sub, 1, mean), apply(proj.sub, 1, var))

curr.sub <- log10(curr.sub+1)
proj.sub <- log10(proj.sub+1)

m <- t(curr.sub)
m <- m - apply(m, 1, mean)
range(apply(m, 1, mean))
m <- m / apply(m, 1, sd)
range(apply(m, 1, var))
foo <- pca$v[,1:nPCs]; rownames(foo) <- ods.genes
pca.curr <- m[,ods.genes2] %*% foo[ods.genes2,]
rownames(pca.curr) <- colnames(curr.sub)
head(pca.curr)

m <- t(proj.sub)
m <- m - apply(m, 1, mean)
range(apply(m, 1, mean))
m <- m / apply(m, 1, sd)
range(apply(m, 1, var))
foo <- pca$v[,1:nPCs]; rownames(foo) <- ods.genes
pca.proj <- m %*% foo[ods.genes2,]
rownames(pca.proj) <- colnames(proj.sub)
head(pca.proj)

colnames(pca.curr) <- colnames(pca.proj) <- paste0('PC', 1:ncol(pca.curr))

## double check zero centered
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

set.seed(0)
gsim.gg = graphViz(t(pca.curr), t(pca.proj), 10,
                   cell.colors=cell.cols,
                   similarity_threshold=0,
                   distance_weight = 1,
                   weighted=FALSE,
                   plot = FALSE,
                   return_graph = TRUE)
par(mfrow=c(1,1))
plot(gsim.gg$fdg_coords, main = "FDG: vertex coordinates", col=cell.cols, pch=16)
#plot(emb_umap, main = "Umap", col=cell.cols, pch=16)

par(mfrow=c(1,1))
set.seed(0)
g <- gsim.gg$graph
V(g)$label = NA
V(g)$size = 2
E(g)$arrow.size = 1
E(g)$color = rgb(0,0,0,0.01)
plot(g)

