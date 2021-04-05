
<a href="https://jef.works/veloviz/"><img src="https://github.com/JEFworks/veloviz/blob/package/docs/img/logo_final.png" width="200"/></a>

<!-- badges: start -->
[![R build status](https://github.com/JEFworks/veloviz/workflows/R-CMD-check/badge.svg)](https://github.com/JEFworks/veloviz/actions)
<!-- badges: end -->


`VeloViz` creates an RNA-velocity-informed 2D embedding for single cell transcriptomics data.

![](https://github.com/JEFworks/veloviz/blob/package/docs/img/readme_schematic.png)

The overall approach is detailed in the [preprint](https://www.biorxiv.org/content/10.1101/2021.01.28.425293v1).

## Installation

To install `VeloViz`, we recommend using `devtools`:

``` r
require(devtools)
devtools::install_github('JEFworks-Lab/veloviz')
```

## Example

Below is a short example showing how to create a VeloViz embedding on sc-RNAseq data.   

``` r
library(veloviz)
# load built in scRNA-seq data
clusters = pancreas$clusters # cell type annotations
pcs = pancreas$pcs # PCs used to make other embeddings (UMAP,tSNE..)
vel = pancreas$vel # RNA velocity

curr = vel$current # current transcriptional state
proj = vel$projected # projected transcriptional state

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
![](https://github.com/JEFworks/veloviz/blob/package/docs/img/readme_example.png)

## Tutorials
[scRNA-seq data preprocessing and visualization using VeloViz](https://github.com/JEFworks/veloviz/blob/package/docs/pancreas.md)  
[MERFISH cell cycle visualization using VeloViz](https://github.com/JEFworks/veloviz/blob/package/docs/merfish.md)  
[Understanding VeloViz parameters](https://github.com/JEFworks/veloviz/blob/package/docs/simulation.md) \
[Visualizing the VeloViz graph using UMAP](https://github.com/JEFworks/veloviz/blob/package/docs/umap.md)