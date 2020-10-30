setwd('graphViz')
library(igraph)
library(matie)
library(RANN)
library(Rcpp)
source("projectedNeighbors_weightedCD.R")

veloviz <- function(mat, curr, proj,
                    cell.cols = NA,
                    center = TRUE,
                    scale = FALSE,
                    alpha = 0.05,
                    nPCs = 10,
                    k = 10,
                    similarity.threshold = 0,
                    distance.weight = 1,
                    weighted = FALSE,
                    cluster.method = igraph::cluster_louvain,
                    layout.method = igraph::layout_with_fr,
                    plot = TRUE,
                    verbose = TRUE,
                    seed = 0
) {

  if(verbose) {
    print('normalizing variance')
  }
  ## identify overdispsersed genes
  ## and establish gene scaling factor
  ## from main full matrix
  matnorm.info <- normalizeVariance(mat, details=TRUE,
                                    plot=FALSE, alpha = alpha, verbose=verbose)
  ods.genes <- rownames(matnorm.info$mat)[matnorm.info$ods]
  ## limit to velocity informative genes
  ods.genes <- intersect(ods.genes, rownames(curr))
  ## scale gene variance in current and future expression
  curr.sub <- curr[ods.genes,]
  proj.sub <- proj[ods.genes,]
  curr.sub <- curr.sub * matnorm.info$df[ods.genes,]$gsf
  proj.sub <- proj.sub * matnorm.info$df[ods.genes,]$gsf
  if(plot) {
    par(mfrow=c(3,2))
    x = log10(apply(mat, 1, mean))
    y = log10(apply(mat, 1, var))
    smoothScatter(x, y, main='Reference (before variance normalization)',
                  xlab = 'log10(mean)', ylab='log10(var)')
    points(x[ods.genes], y[ods.genes], col='red', pch='.')

    mat.sub <- mat * matnorm.info$df$gsf
    x = log10(apply(mat.sub, 1, mean))
    y = log10(apply(mat.sub, 1, var))
    smoothScatter(x, y, main='Reference (after variance normalization)',
                  xlab = 'log10(mean)', ylab='log10(var)')
    points(x[ods.genes], y[ods.genes], col='red', pch='.')

    x = log10(apply(curr, 1, mean))
    y = log10(apply(curr, 1, var))
    smoothScatter(x, y, main='Current (before variance normalization)',
                  xlab = 'log10(mean)', ylab='log10(var)')
    points(x[ods.genes], y[ods.genes], col='red', pch='.')

    x = log10(apply(curr.sub, 1, mean))
    y = log10(apply(curr.sub, 1, var))
    smoothScatter(x, y, main='Current (after variance normalization)',
                  xlab = 'log10(mean)', ylab='log10(var)')

    x = log10(apply(proj, 1, mean))
    y = log10(apply(proj, 1, var))
    smoothScatter(x, y, main='Current (before variance normalization)',
                   xlab = 'log10(mean)', ylab='log10(var)')
    points(x[ods.genes], y[ods.genes], col='red', pch='.')

    x = log10(apply(proj.sub, 1, mean))
    y = log10(apply(proj.sub, 1, var))
    smoothScatter(x, y, main='Current (after variance normalization)',
                   xlab = 'log10(mean)', ylab='log10(var)')
  }
  ## log transform
  curr.sub <- log10(curr.sub+1)
  proj.sub <- log10(proj.sub+1)

  if(verbose) {
    print('dimensionality reduction by PCA')
  }
  ## establish PCs from overdispersed genes
  ## on main full matrix(?)
  m <- log10(matnorm.info$mat[ods.genes,]+1)
  ## on current matrix(?)
  #m <- curr.sub
  pca <- RSpectra::svds(A = t(m), k=min(50, nPCs+10),
              opts = list(center = center, scale = scale,
                          maxitr = 2000, tol = 1e-10))

  if(plot) {
    par(mfrow=c(1,1))
    plot(pca$d)
    abline(v=nPCs, col='red')
  }

  ## project current cells onto PCs
  m <- t(curr.sub)
  if(center) {
    m <- m - apply(m, 1, mean)
  }
  if(scale) {
    m <- m / apply(m, 1, sd)
  }
  pca.curr <- m %*% pca$v[,1:nPCs]
  rownames(pca.curr) <- colnames(curr.sub)

  ## project future onto PCs
  m <- t(proj.sub)
  if(center) {
    m <- m - apply(m, 1, mean)
  }
  if(scale) {
    m <- m / apply(m, 1, sd)
  }
  pca.proj <- m %*% pca$v[,1:nPCs]
  rownames(pca.proj) <- colnames(proj.sub)
  colnames(pca.curr) <- colnames(pca.proj) <- paste0('PC', 1:ncol(pca.curr))

  if(verbose) {
    print('generating velocity informed embedding')
  }
  set.seed(seed)
  gsim.gg <- graphViz(t(pca.curr), t(pca.proj), k,
                      cell.colors=NA,
                      similarity_threshold=similarity.threshold,
                      distance_weight = distance.weight,
                      weighted = weighted,
                      plot = FALSE,
                      return_graph = TRUE)

  if(verbose) {
    ## other useful statistics to evaluate graph stability?
    print(paste0('graph transivity: ', igraph::transitivity(gsim.gg$graph)))
  }

  if(is.na(cell.cols)){
    print('identifying clusters')
    g <- gsim.gg$graph
    g <- igraph::simplify(g)
    g <- igraph::as.undirected(g)
    km <- cluster.method(g)
    com <- km$membership
    names(com) <- km$names
    com <- factor(com)
    col = rainbow(length(levels(com)),s = 0.8, v = 0.8)
    cell.cols = col[com] #color according to cluster
    names(cell.cols) = names(com)
  }

  if(plot) {
    par(mfrow=c(1,1))
    #plot(gsim.gg$fdg_coords, main = "FDG: vertex coordinates", col=cell.cols, pch=16)

    set.seed(seed)
    g <- gsim.gg$graph
    V(g)$label = NA
    V(g)$size = 2
    V(g)$color = cell.cols
    V(g)$frame.color = NA
    E(g)$arrow.size = 0.5
    E(g)$color = rgb(0,0,0,0.05)
    coords <- layout.method(g)
    plot(g, layout = coords)
  }

  return(gsim.gg)
}

############ Pancreas example
## try on example
library(reticulate)
library(Matrix)
conda_list()
use_condaenv("r-velocity", required = TRUE)
scv = import("scvelo")

adata = scv$datasets$pancreas()
spliced = as.matrix(t(adata$layers['spliced']))
unspliced = as.matrix(t(adata$layers['unspliced']))
cells = adata$obs_names$values
genes = adata$var_names$values
colnames(spliced) = colnames(unspliced) = cells
rownames(spliced) = rownames(unspliced) = genes
mat <- spliced+unspliced
vi <- rowSums(mat) > 100
table(vi)
mat <- mat[vi,]

clusters = adata$obs$clusters #extract clusters
names(clusters) = adata$obs_names$values
col = rainbow(length(levels(clusters)),s = 0.8, v = 0.8)
cell.cols = col[clusters] #color according to cluster
names(cell.cols) = names(clusters)

library(MUDAN) ## just port over functions later
vel = readRDS('../panc_vel_k30.rds')
#panc_k30 <- veloviz(mat, vel$current, vel$projected, cell.cols)
## use convolved matrix because
## it looks like the current and projected expressions are already convolved
## so the variance scaling doesn't make sense if we build it on a main matrix
## that doesn't match the curr and proj
panc_k30 <- veloviz(vel$conv.emat.norm + vel$conv.nmat.norm,
                    vel$current, vel$projected, cell.cols)
## or should we just use the current matrix?
## do we expect the mean variance dependency made on all genes to be so different?
panc_k30 <- veloviz(vel$current,
                    vel$current, vel$projected, cell.cols)
#vel = readRDS('../panc_vel_k1.rds')
#panc_k1 <- veloviz(vel$conv.emat.norm + vel$conv.nmat.norm, vel$current, vel$projected, cell.cols)
#vel = readRDS('../panc_vel_k100.rds')
#panc_k100 <- veloviz(vel$conv.emat.norm + vel$conv.nmat.norm, vel$current, vel$projected, cell.cols)

######### merfish example
load(file="../merfish_vel_k30.RData")
vel = rvel.cd
## variance looks a bit weird
merfish_k30 <- veloviz(mat = vel$conv.emat.norm + vel$conv.nmat.norm,
                       curr = vel$current, proj = vel$projected, cell.cols=cell.colors, nPCs = 10, k = 50)
## gene scaling factors look more reasonable
merfish_k30 <- veloviz(mat = vel$current,
                    curr = vel$current, proj = vel$projected, cell.cols=cell.colors, nPCs = 10, k = 50)
## is transivity helpful?
merfish_k30 <- veloviz(mat = vel$current,
                       curr = vel$current, proj = vel$projected, cell.cols=cell.colors, nPCs = 3, k = 10)

########## Neuro example
vel = readRDS(file="../neuro_vel_k30.rds")
neuro_k30 <- veloviz(vel$conv.emat.norm + vel$conv.nmat.norm,
                    vel$current, vel$projected)
neuro_k30 <- veloviz(vel$current,
                    vel$current, vel$projected)
