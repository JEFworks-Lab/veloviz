#' Normalizes counts to CPM
#'
#' Normalizes raw counts to counts per million
#'
#' @param counts Read count matrix. The rows correspond to genes, columns
#' correspond to individual cells
#' @param depthScale Depth scaling. Using a million for CPM (default: 1e6)
#' @param verbose Boolean for verbosity setting (default: TRUE)
#'
#' @return a normalized matrix
#'
#' @export
normalizeDepth <- function(
  counts                ,
  depthScale     = 1e+06,
  verbose        = TRUE
){
  if (!class(counts) %in% c("dgCMatrix", "dgTMatrix")) {
    if (verbose) {
      message("Converting to sparse matrix ...")
    }
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  if (verbose) {
    message("Normalizing matrix with ", ncol(counts), " cells and ",
            nrow(counts), " genes")
  }
  counts <- Matrix::t(Matrix::t(counts)/Matrix::colSums(counts))
  counts <- counts * depthScale
  return(counts)
}


#' Identify overdispersed genes by normalizing counts per million (CPM)
#' gene expression variance relative to transcriptome-wide expectations
#' (Modified from SCDE/PAGODA2 code)
#'
#' Normalizes gene expression magnitudes to with respect to its ratio to the
#' transcriptome-wide expectation as determined by local regression on
#' expression magnitude
#'
#' @param cpm Counts per million (CPM) matrix. Rows are genes, columns are cells.
#' @param gam.k Generalized additive model parameter; the dimension of the
#' basis used to represent the smooth term (default: 5)
#' @param alpha Significance threshold for overdispersed genes (default: 0.05)
#' @param max.adjusted.variance Ceiling on maximum variance after normalization
#' to prevent infinites (default: 1e3)
#' @param min.adjusted.variance Floor on minimum variance after normalization
#' (default: 1e-3)
#' @param verbose Boolean for verbosity setting (default: TRUE)
#' @param details Boolean to return detailed data frame or names of genes
#' (default: FALSE)
#'
#' @return A list with two items: (1) an adjusted CPM matrix with the same
#' dimensions as the input and (2) a dataframe with the summary statistics for
#' each gene.
#'
#' @export
normalizeVariance <- function(
  cpm,
  gam.k                 = 5,
  alpha                 = 0.05,
  max.adjusted.variance = 5,
  min.adjusted.variance = 1,
  verbose               = TRUE,
  plot                  = FALSE,
  details               = FALSE
) {

  if (!class(cpm) %in% c("dgCMatrix", "dgTMatrix")) {
    if (verbose) {
      message("Converting to sparse matrix ...")
    }
    cpm <- Matrix::Matrix(cpm, sparse = TRUE)
  }

  n.cells <- ncol(cpm)
  n.obs <- nrow(cpm)

  ## faster to use moments to calculate variance
  ## for sparse matrix
  cpm_center <- Matrix::rowMeans(cpm)
  sumx      <- Matrix::rowSums(cpm)
  sumxx     <- Matrix::rowSums(cpm^2)
  cpm_var   <- sumxx - 2 * sumx * cpm_center + ncol(cpm) * cpm_center ^ 2

  df <- data.frame(
    log_mean     = log(cpm_center),
    log_variance = log(cpm_var)
  )
  rownames(df) <- rownames(cpm)

  vi <- which(is.finite(df$log_variance))

  # Too few genes
  if (length(vi) < gam.k * 1.5) {
    gam.k <- 1
  }

  if (gam.k < 2) {
    if(verbose) {
      message("Using linear modeling...")
    }
    m <- lm(log_variance ~ log_mean, data = df[vi,])
  } else {
    if(verbose) {
      message(paste0("Using general additive modeling with k = ", gam.k, "..."))
    }
    s <- mgcv::s
    fm <- as.formula(sprintf("log_variance ~ s(log_mean, k = %s)", gam.k))
    m <- mgcv::gam(fm, data = df[vi,])
  }

  df$res <- -Inf
  df$res[vi] <- resid(m, type = "response")

  df$log_p  <- as.numeric(pf(
    q          = exp(df$res),
    df1        = n.obs,
    df2        = n.obs,
    lower.tail = FALSE,
    log.p      = TRUE
  ))
  df$log_p_adjusted <- bh.adjust(df$log_p, log=TRUE)

  df$qv  <- as.numeric(
    qchisq(
      p          = df$log_p,
      df         = n.cells - 1,
      lower.tail = FALSE,
      log.p      = TRUE
    ) / n.cells
  )

  df$over_disp <- df$log_p_adjusted < log(alpha)
  if(verbose) {
    message(paste0("Identifed ", sum(df$over_disp), " overdispersed genes using
    adjusted p-value threshold alpha = ", alpha))
  }

  clamp <- function(x, min, max) pmax(min, pmin(max, x))

  #df$scale_factor <- sqrt(
  #  clamp(df$qv, min.adjusted.variance, max.adjusted.variance) /
  #    exp(df$log_variance)
  #)
  #df$scale_factor <- clamp(df$res, min.adjusted.variance, max.adjusted.variance)
  df$scale_factor <- clamp(df$qv, min.adjusted.variance, max.adjusted.variance)
  df$scale_factor[!is.finite(df$scale_factor)] <- 0

  if(plot) {
    par(mfrow=c(1,2))
    plot(df$log_mean, df$log_var,
         col = ifelse(df$over_disp, "red", "black"), pch=".")
    plot(df$log_mean, df$res,
         col = ifelse(df$over_disp, "red", "black"), pch=".")
    par(mfrow=c(1,1))
    plot(df$res, df$scale_factor,
         col = ifelse(df$over_disp, "red", "black"), pch=".")
  }

  if(details) {
    return(df)
  } else {
    return(rownames(df)[df$over_disp])
  }
}
## BH P-value adjustment with a log option
bh.adjust <- function(x, log = FALSE) {
  nai <- which(!is.na(x))
  ox <- x
  x <- x[nai]
  id <- order(x, decreasing = FALSE)
  if(log) {
    q <- x[id] + log(length(x)/seq_along(x))
  } else {
    q <- x[id]*length(x)/seq_along(x)
  }
  a <- rev(cummin(rev(q)))[order(id)]
  ox[nai] <- a
  ox
}

## PCA
reduceDimensions <- function(cpm,
                             center=TRUE,
                             scale=FALSE,
                             use.ods.genes=TRUE,
                             max.ods.genes=1000,
                             alpha=0.05,
                             nPCs=50,
                             verbose=TRUE,
                             plot=FALSE,
                             details=FALSE)
{
  ods.genes = normalizeVariance(cpm, alpha=alpha, verbose=verbose, plot=plot, details = TRUE)
  scale.factor = exp(ods.genes$scale_factor) ## use residual variance q-value as scale factor
  names(scale.factor) <- rownames(cpm)

  if(use.ods.genes) {
    if(verbose) {
      print('Normalizing variance...')
    }
    if(sum(ods.genes$over_disp) > max.ods.genes) {
      if(verbose) {
        print(paste0('Limiting to top ', max.ods.genes, ' overdispersed genes...'))
      }
      best.genes <- names(sort(scale.factor, decreasing=TRUE)[1:max.ods.genes])
      cpm = cpm[best.genes,]
      scale.factor = scale.factor[best.genes]
    } else {
      cpm = cpm[ods.genes$over_disp,]
      scale.factor = scale.factor[ods.genes$over_disp]
    }
  }

  ## establish PCs from overdispersed genes
  m <- cpm
  ## mean
  rmean <- Matrix::rowMeans(m)
  sumx     <- Matrix::rowSums(m)
  sumxx    <- Matrix::rowSums(m^2)
  ## sd
  rsd <- sqrt(sumxx - 2 * sumx * rmean + ncol(m) * rmean ^ 2)
  if(center) {
    m_center <- rmean
  } else {
    m_center <- FALSE
  }
  if(scale) {
    ## regular scale to var = 1
    m_scale <- rsd
  } else {
    ## scale to residual variance
    m_scale <- rsd/scale.factor
  }
  pca <- RSpectra::svds(A = Matrix::t(m),
                        k=min(50, nPCs+10),
                        opts = list(
                          center = m_center,
                          scale = m_scale,
                          maxitr = 2000,
                          tol = 1e-10))

  if(plot) {
    plot(pca$d, type="l")
    abline(v=nPCs, col='red')
  }

  pcs <- pca$u[,1:nPCs]
  rownames(pcs) <- colnames(cpm)
  colnames(pcs) <- paste0('PC', 1:nPCs)

  if(details) {
    return(pca)
  } else {
    return(as.matrix(pcs))
  }
}

