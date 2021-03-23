#' Creates VeloViz graph and FDG layout from PC scores of current and projected transcriptional states.  
#'
#' @param curr Genes (rows) x cells (columns) matrix of observed current transcriptional state 
#' @param proj Genes (rows) x cells (columns) matrix of predicted future transcriptional state 
#' @param normalize.depth logical to normalize raw counts to counts per million, default = TRUE
#' @param depth Depth scaling, default = 1e6 for counts per million (CPM)
#' @param use.ods.genes Use only overdispersed genes to create VeloViz graph, default = TRUE
#' @param max.ods.genes number of most highly expressed overdispersed genes to use to create VeloViz graph, default = 2000
#' @param alpha Significance threshold for overdispersed genes, default = 0.05
#' @param pca logical to use PC scores to create VeloViz graph, default = TRUE. FALSE = use gene expression to create VeloViz graph
#' @param center logical to mean center gene expression before PCA, default = TRUE
#' @param scale logical to scale gene expression variance before PCA, default = TRUE
#' @param nPCs number of principal components to use to create VeloViz graph, default = 10
#' @param k Number of nearest neighbors to assign each cell
#' @param similarity.threshold similarity threshold below which to remove edges, default = -1 i.e. no edges removed
#' @param distance.weight Weight of distance component of composite distance, default = 1
#' @param distance.threshold quantile threshold for distance component above which to remove edges, default = 1 i.e. no edges removed
#' @param weighted logical indicating whether to compute VeloViz edges based on composite distance, default = TRUE. FALSE = all edges are of equal weight
#' @param remove.unconnected logical indicating whether to remove cells with no edges in the VeloViz graph from the output embedding, default = TRUE (removed)
#' @param verbose logical for verbosity setting, default = FALSE
#' @param details logical to return detailed data frame or names of genes, default = FALSE
#' @param seed seed to supply FDG function for reproducible layout
#' 
#' @return `graph` igraph object of VeloViz graph
#' @return `fdg_coords` cells (rows) x 2 coordinates of force-directed layout of VeloViz graph
#' @return `projectedNeighbors` output of `projectedNeighbors`
#' 
#' @examples 
#' vel <- pancreas$vel
#' curr <- vel$current
#' proj <- vel$projected
#' 
#' buildVeloviz(curr = curr, proj = proj, normalize.depth = TRUE, 
#' use.ods.genes = TRUE, alpha = 0.05, pca = TRUE, nPCs = 20, center = TRUE, 
#' scale = TRUE, k = 5, similarity.threshold = 0.25, distance.weight = 1,
#' distance.threshold = 0.5, weighted = TRUE, seed = 0, verbose = FALSE)
#' 
#' @seealso \code{\link{projectedNeighbors}}
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
                         distance.threshold = 1,
                         weighted = TRUE,
                         remove.unconnected = TRUE,
                         verbose = FALSE,
                         details = FALSE,
                         seed = 0
) {

  if (!class(curr)[1] %in% c("dgCMatrix", "dgTMatrix")) {
    if (verbose) {
      message("Converting to sparse matrix ...")
    }
    curr <- Matrix::Matrix(curr, sparse = TRUE)
  }
  if (!class(proj)[1] %in% c("dgCMatrix", "dgTMatrix")) {
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
    
    ## check if curr or proj have negative values 
    if(sum(curr<0)>0){
      stop('curr contains negative values, cannot log normalize for PCA')
    }
    
    if(sum(proj<0)>0){
      stop('proj contains negative values, cannot log normalize for PCA')
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
    pca.curr = t(curr)
    pca.proj = t(proj)
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
                      distance_threshold = distance.threshold,
                      weighted = weighted,
                      remove_unconnected = remove.unconnected,
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
#' @param vig output of buildVeloviz
#' @param layout.method igraph method to use for generating 2D graph representation, default = igraph::layout_with_fr
#' @param clusters cluster annotations for cells in data 
#' @param cluster.method igraph method to use for clustering if clusters are not provided, default = igraph::cluster_louvain
#' @param col colors to use for plotting
#' @param alpha transparency for plotting graph nodes
#' @param verbose logical for verbosity setting, default = FALSE
#' @param seed seed to supply FDG function for reproducible layout
#'
#' @return cells (rows) x 2 coordinates of force-directed layout of VeloViz graph
#'
#' @examples 
#' vel <- pancreas$vel
#' curr <- vel$current
#' proj <- vel$projected
#' 
#' vv <- buildVeloviz(curr = curr, proj = proj, normalize.depth = TRUE, 
#' use.ods.genes = TRUE, alpha = 0.05, pca = TRUE, nPCs = 20, center = TRUE, 
#' scale = TRUE, k = 5, similarity.threshold = 0.25, distance.weight = 1,
#' distance.threshold = 0.5, weighted = TRUE, seed = 0, verbose = FALSE)
#' 
#' plotVeloviz(vv)
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
