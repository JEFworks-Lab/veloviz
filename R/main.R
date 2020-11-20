#' Main function
#'
#' @export
#'
buildVeloviz <- function(curr, proj,
                         normalize.depth = TRUE,
                         depth = 1e6,
                         use.ods.genes = TRUE,
                         max.ods.genes = 2000,
                         alpha = 0.05,
                         pca = TRUE,
                         center = TRUE,
                         scale = TRUE,
                         nPCs = 10,
                         k = 10,
                         similarity.threshold = 0,
                         distance.weight = 1,
                         weighted = FALSE,
                         verbose = FALSE,
                         details = FALSE,
                         seed = 0
) {

  if (!class(curr) %in% c("dgCMatrix", "dgTMatrix")) {
    if (verbose) {
      message("Converting to sparse matrix ...")
    }
    curr <- Matrix::Matrix(curr, sparse = TRUE)
  }
  if (!class(proj) %in% c("dgCMatrix", "dgTMatrix")) {
    if (verbose) {
      message("Converting to sparse matrix ...")
    }
    proj <- Matrix::Matrix(proj, sparse = TRUE)
  }

  if(normalize.depth) {
    if(verbose) {
      message('Normalizing depth...')
    }
    curr = normalizeDepth(curr, depthScale = depth, verbose=verbose)
    proj = normalizeDepth(proj, depthScale = depth, verbose=verbose)
  }

  if(use.ods.genes) {
    if(verbose) {
      message('Identifying overdispersed genes...')
    }
    matnorm.info = normalizeVariance(curr, alpha=alpha, verbose=verbose, details = TRUE)
    curr = matnorm.info$matnorm
    scale.factor = matnorm.info$df$scale_factor
    names(scale.factor) <- rownames(matnorm.info$df)

    ## scale variance
    m <- proj
    rmean <- Matrix::rowMeans(m)
    sumx     <- Matrix::rowSums(m)
    sumxx    <- Matrix::rowSums(m^2)
    rsd <- sqrt((sumxx - 2 * sumx * rmean + ncol(m) * rmean ^ 2) / (ncol(m)-1))
    ## use same scale factor as curr
    proj <- proj / rsd * scale.factor[names(rsd)]
    proj <- proj[rownames(curr),]

    if(nrow(curr) > max.ods.genes) {
      if(verbose) {
        message(paste0('Limiting to top ', max.ods.genes, ' highly expressed genes...'))
      }
      rmean <- Matrix::rowMeans(curr)
      best.genes <- names(sort(rmean, decreasing=TRUE)[1:max.ods.genes])
      curr = curr[best.genes,]
      proj = proj[best.genes,]
    }
  }

  if(pca) {
    if(verbose) {
      message('Performing dimensionality reduction by PCA...')
    }
    ## establish PCs from overdispersed genes
    m <- log10(curr+1)
    ## mean
    rmean <- Matrix::rowMeans(m)
    sumx     <- Matrix::rowSums(m)
    sumxx    <- Matrix::rowSums(m^2)
    ## sd
    rsd <- sqrt((sumxx - 2 * sumx * rmean + ncol(m) * rmean ^ 2) / (ncol(m)-1))
    if(center) {
      if(verbose) {
        message('Centering...')
      }
      m <- m - rmean
    }
    if(scale) {
      if(verbose) {
        message('Using unit variance...')
      }
      ## regular scale to var = 1
     m <- m/rsd
    }
    pca <- RSpectra::svds(A = Matrix::t(m),
                          k=min(50, nPCs+10),
                          opts = list(
                            center = FALSE, ## already done
                            scale = FALSE, ## already done
                            maxitr = 2000,
                            tol = 1e-10))

    ## project current cells onto PCs
    if(verbose) {
      message('Projecting current cells onto PCs...')
    }
    m <- log10(curr+1)
    if(center) {
      m <- m - rmean
    }
    if(scale) {
      m <- m / rsd
    }
    pca.curr <- Matrix::t(m) %*% pca$v[,1:nPCs]

    ## project future onto PCs
    if(verbose) {
      message('Projecting future cells onto PCs...')
    }
    m <- log10(proj+1)
    rmean <- Matrix::rowMeans(m)
    sumx     <- Matrix::rowSums(m)
    sumxx    <- Matrix::rowSums(m^2)
    ## sd
    rsd <- sqrt((sumxx - 2 * sumx * rmean + ncol(m) * rmean ^ 2) / (ncol(m)-1))
    if(center) {
      m <- m - rmean
    }
    if(scale) {
      m <- m / rsd
    }
    pca.proj <- Matrix::t(m) %*% pca$v[,1:nPCs]

    colnames(pca.curr) <- colnames(pca.proj) <- paste0('PC', 1:ncol(pca.curr))
  } else {
    pca.curr = curr
    pca.proj = proj
  }

  if(verbose) {
    message('Generating velocity informed embedding...')
  }

  ## dense matrix by now
  pca.curr <- as.matrix(pca.curr)
  pca.proj <- as.matrix(pca.proj)

  set.seed(seed)
  ## velocity informed graph
  vig <- graphViz(t(pca.curr), t(pca.proj), k,
                      cell.colors=NA,
                      similarity_threshold=similarity.threshold,
                      distance_weight = distance.weight,
                      weighted = weighted,
                      plot = FALSE,
                      return_graph = TRUE)

  if(details) {
    return(list(
      pca.curr = pca.curr,
      pca.proj = pca.proj,
      vig = vig
    ))
  } else {
    return(vig)
  }
}

#' Plot function
#'
#' @export
#'
plotVeloviz <- function(
  vig,
  layout.method = igraph::layout_with_fr,
  clusters = NA,
  cluster.method = igraph::cluster_louvain,
  col = NA,
  alpha = 0.05,
  verbose = TRUE,
  seed = 0
) {

  if(!is.na(clusters) & is.na(col)) {
    if(verbose) {
      message('Using provided clusters...')
    }
    com <- as.factor(clusters)
    col = rainbow(length(levels(com)), s = 0.8, v = 0.8)
    cell.cols <- col[clusters]
  } else if(is.na(clusters) & !is.na(col)) {
    if(verbose) {
      message('Using provided colors...')
    }
    cell.cols <- col
  } else if(is.na(clusters & is.na(col))){
    if(verbose) {
      message('Identifying clusters...')
    }
    g <- vig$graph
    g <- igraph::simplify(g)
    g <- igraph::as.undirected(g)
    km <- cluster.method(g)
    com <- km$membership
    names(com) <- km$names
    com <- factor(com)
    col = rainbow(length(levels(com)), s = 0.8, v = 0.8)
    cell.cols = col[com] #color according to cluster
    names(cell.cols) = names(com)
  }

  set.seed(seed)
  g <- vig$graph
  igraph::V(g)$label = NA
  igraph::V(g)$size = 2
  igraph::V(g)$color = cell.cols
  igraph::V(g)$frame.color = NA
  igraph::E(g)$arrow.size = 0.5
  igraph::E(g)$color = rgb(0,0,0,0.05)
  coords <- layout.method(g)
  rownames(coords) <- rownames(vig$fdg_coords)
  plot(g, layout = coords)

  return(coords)
}
