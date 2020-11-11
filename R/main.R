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
