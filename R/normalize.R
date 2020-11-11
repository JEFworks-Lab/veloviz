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
#' @importFrom Matrix Matrix colSums t
normalizeCounts <- function(
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
#' @importFrom mgcv s
#' @export
getOverdispersedGenes <- function(
  cpm,
  gam.k                 = 5,
  alpha                 = 0.05,
  max.adjusted.variance = 1e3,
  min.adjusted.variance = 1e-3,
  verbose               = TRUE,
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
      message(paste0("Using general additive modeling with k = ", k, "..."))
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

  df$scale_factor <- sqrt(
    clamp(df$qv, min.adjusted.variance, max.adjusted.variance) /
      exp(df$log_variance)
  )
  df$scale_factor[!is.finite(df$scale_factor)] <- 0

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

