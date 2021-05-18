#' Pancreas scRNA-seq data
#' 
#' Pancreatic endocrinogenesis scRNA-seq from Bastidas-Ponce et. al., 
#' Development 2019 accessed via scVelo package and 
#' subsampled to 739 cells. 
#' 
#' 
#' 
#' @format list of 4 objects:
#' \describe{
#' \item{spliced}{matrix, 7192 genes x 739 cells of spliced counts}
#' \item{unspliced}{matrix, 7192 genes x 739 cells of unspliced counts}
#' \item{pcs}{matrix, 739 x 50 cell scores in 50 PCs}
#' \item{clusters}{factor of cell type annotations from scVelo}
#' 
#' }
#'
#' @source \url{https://dev.biologists.org/content/146/12/dev173849.long}
"pancreas"

#' MERFISH velocity subset
#' 
#' output of velocyto.R::gene.relative.velocity.estimates for 40 cell subset of MERFISH data. Used to run examples
#' 
#' 
#' 
#' @format list of 1:
#' \describe{
#' \item{vel}{velocity output containing current observed ("current") and predicted future ("projected") estimates}
#' }
#'
#' @source \url{https://www.pnas.org/content/116/39/19490}
"vel"
