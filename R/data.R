#' Pancreas scRNA-seq data
#' 
#' Pancreatic endocrinogenesis scRNA-seq from Bastidas-Ponce et. al., 
#' Development 2019 accessed via scVelo package and 
#' subsampled to 739 cells. 
#' 
#' 
#' 
#' @format list of 6 objects:
#' \describe{
#' \item{spliced}{matrix, 7192 genes x 739 cells of spliced counts}
#' \item{unspliced}{matrix, 7192 genes x 739 cells of unspliced counts}
#' \item{clusters}{factor of cell type annotations from scVelo}
#' \item{pcs}{matrix, 739 x 50 cell scores in 60 PCs}
#' \item{cell.dist}{dist, pairwise cell distances in PC space used to compute velocity}
#' \item{vel}{list, output from running velocyto using spliced, unspliced, and cell.dist}
#' }
#'
#' @source \url{https://dev.biologists.org/content/146/12/dev173849.long}
"pancreas"

#' Pancreas scRNA-seq data missing intermediates
#' 
#' Pancreatic endocrinogenesis scRNA-seq from Bastidas-Ponce et. al., 
#' Development 2019 accessed via scVelo package with missing Ngn3 high EP cells
#' and subsampled to 660 cells. 
#'
#' @format list of 6 objects:
#' \describe{
#' \item{spliced}{matrix, 7192 genes x 739 cells of spliced counts}
#' \item{unspliced}{matrix, 7192 genes x 739 cells of unspliced counts}
#' \item{clusters}{factor of cell type annotations from scVelo}
#' \item{pcs}{matrix, 739 x 50 cell scores in 60 PCs}
#' \item{cell.dist}{dist, pairwise cell distances in PC space used to compute velocity}
#' \item{vel}{list, output from running velocyto using spliced, unspliced, and cell.dist}
#' }
#'
#' @source \url{https://dev.biologists.org/content/146/12/dev173849.long}
"pancreasWithGap"

#' MERFISH scRNA-seq data
#' 
#' MERFISH sequencing of U20S cells in culture, subsampled to 645 cells. 
#' 
#' 
#' 
#' @format list of 6 objects:
#' \describe{
#' \item{nuc}{matrix, 9050 genes x 645 cells of nuclear counts}
#' \item{cyto}{matrix, 9050 genes x 645 cells of cytoplasmic counts}
#' \item{col}{factor of cell colors corresponding to position in cell cycle}
#' \item{pcs}{matrix, 645 x 50 cell scores in 50 PCs}
#' \item{cell.dist}{dist, pairwise cell distances in PC space used to compute velocity}
#' \item{vel}{list, output from running velocyto using spliced, unspliced, and cell.dist}
#' }
#'
#' @source \url{https://www.pnas.org/content/116/39/19490}
"MERFISH"
