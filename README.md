
<a href="https://jef.works/veloviz/"><img src="https://github.com/JEFworks-Lab/veloviz/blob/package_extras/docs/img/logo_final.png" width="200"/></a>

<!-- badges: start -->
[![R build status](https://github.com/JEFworks/veloviz/workflows/R-CMD-check/badge.svg)](https://github.com/JEFworks/veloviz/actions)
<!-- badges: end -->


`VeloViz` creates an RNA-velocity-informed 2D embedding for single cell transcriptomics data.

![](https://github.com/JEFworks-Lab/veloviz/blob/package_extras/docs/img/readme_schematic.png)

The overall approach is detailed in the [preprint](https://www.biorxiv.org/content/10.1101/2021.01.28.425293v2).

## Installation

To install `VeloViz`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks-Lab/veloviz')
```

VeloViz is also available on [Bioconductor](https://bioconductor.org/packages/veloviz):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("veloviz")

```


## Example

Below is a short example showing how to create a VeloViz embedding on sc-RNAseq data. R objects containing the preprocessed data and velocity models used in this example can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.4632471).   

``` r
library(veloviz)
library(velocyto.R)
# load built in pancreas scRNA-seq data
data(pancreas)
spliced = pancreas$spliced
unspliced = pancreas$unspliced
clusters = pancreas$clusters # cell type annotations
pcs = pancreas$pcs # PCs used to make other embeddings (UMAP,tSNE..)

#choose colors based on clusters for plotting later
cell.cols = rainbow(8)[as.numeric(clusters)]
names(cell.cols) = names(clusters)

# compute velocity using velocyto.R
#cell distance in PC space
cell.dist = as.dist(1-cor(t(pcs))) # cell distance in PC space

vel = gene.relative.velocity.estimates(spliced,
                                       unspliced,
                                       kCells = 30,
                                       cell.dist = cell.dist,
                                       fit.quantile = 0.1)

curr = vel$current
proj = vel$projected

# build VeloViz graph
veloviz = buildVeloviz(
  curr = curr, proj = proj,
  normalize.depth = TRUE,
  use.ods.genes = TRUE,
  alpha = 0.05,
  pca = TRUE,
  nPCs = 20,
  center = TRUE,
  scale = TRUE,
  k = 5,
  similarity.threshold = 0.25,
  distance.weight = 1,
  distance.threshold = 0.5,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

# extract VeloViz embedding
emb.veloviz = veloviz$fdg_coords

# compare to other embeddings

par(mfrow = c(2,2))
#PCA
emb.pca = pcs[,1:2]
plotEmbedding(emb.pca, groups=pancreas$clusters, main='PCA')

#tSNE
set.seed(0)
emb.tsne = Rtsne::Rtsne(pcs, perplexity=30)$Y
rownames(emb.tsne) = rownames(pcs)
plotEmbedding(emb.tsne, groups=pancreas$clusters, main='tSNE',
              xlab = "t-SNE X", ylab = "t-SNE Y")

#UMAP
set.seed(0)
emb.umap = uwot::umap(pcs, min_dist = 0.5)
rownames(emb.umap) <- rownames(pcs)
plotEmbedding(emb.umap, groups=pancreas$clusters, main='UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")

#VeloViz
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz')

```
![](https://github.com/JEFworks-Lab/veloviz/blob/package_extras/docs/img/readme_example.png)

## Tutorials
[scRNA-seq data preprocessing and visualization using VeloViz](https://github.com/JEFworks-Lab/veloviz/blob/package_extras/docs/pancreas.md)  
[MERFISH cell cycle visualization using VeloViz](https://github.com/JEFworks-Lab/veloviz/blob/package_extras/docs/merfish.md)  
[Understanding VeloViz parameters](https://github.com/JEFworks-Lab/veloviz/blob/package_extras/docs/simulation.md) \
[Visualizing the VeloViz graph using UMAP](https://github.com/JEFworks-Lab/veloviz/blob/package_extras/docs/umap.md) \
[VeloViz with dynamic velocity estimates from scVelo](https://github.com/JEFworks-Lab/veloviz/blob/package_extras/docs/scVeloVignette.md)
